function ret = sfb_ctrl(action, task, num_bandits, varargin)
    default_draw = false;
    default_draw_sfb = false;
    default_epsilon = 0.1;
    default_sparse_decay = 0.99;
    default_sparse_history = 1e3;
    default_sparse_lambda = 1e-3;
    default_max_eigens = 6;
    default_chipwidth = 7;
    default_chipheight = 6;
    default_sparse_file = [];
    default_sparse_save_interval = 100;
    default_sparse_save = false;

    parser = inputParser;
    parser.addRequired('action', @(x)isa(x, 'function_handle'));
    parser.addRequired('task');
    parser.addRequired('num_bandits', @isnumeric);
    parser.addParameter('epsilon', default_epsilon, @isnumeric);
    parser.addParameter('sparse_decay', default_sparse_decay, @isnumeric);
    parser.addParameter('sparse_lambda', default_sparse_lambda, @isnumeric);
    parser.addParameter('sparse_history', default_sparse_history, @isnumeric);
    parser.addParameter('sparse_file', default_sparse_file);
    parser.addParameter('sparse_save', default_sparse_save);
    parser.addParameter('sparse_save_interval', default_sparse_save_interval);
    parser.addParameter('max_eigens', default_max_eigens, @isnumeric);
    parser.addParameter('chipwidth', default_chipwidth);
    parser.addParameter('chipheight', default_chipheight);
    parser.addParameter('draw', default_draw);
    parser.addParameter('draw_sfb', default_draw_sfb);
    parser.KeepUnmatched = true;
    parse(parser, action, task, num_bandits, varargin{:});
    param = parser.Results;
    logging.message('%s\n%s', mfilename, third_party.struct2str(parser.Results));

    if param.draw_sfb
        figure;
        perrow = ceil(sqrt(num_bandits));
        [imx, imy] = ndgrid(0:0.05:1, 0:0.05:1);
        imobjs = cell(1, num_bandits);
        for ind = 1:num_bandits
            x = mod(ind - 1, perrow);
            y = floor((ind - 1) / perrow);
            axes('position', [x / perrow 1 - (y + 1) / perrow 1 / perrow 1 / perrow]);
            imobjs{ind} = image(zeros(size(imx)));
            set(gca, 'visible', 'off');
            colormap hot;
        end
    end

    ret = @step;
    if ~isempty(param.sparse_file) && exist(param.sparse_file, 'file')
        logging.message('Loading sparse fourier model from file ''%s''', param.sparse_file);
        D = load(param.sparse_file);
        p_hist = D.p_hist;
        dp_hist = D.dp_hist;
        weights = D.weights;
    else
        logging.message('Starting modeling from scratch');
        p_hist = cell(1, num_bandits);
        dp_hist = cell(1, num_bandits);
        weights = cell(1, num_bandits);
    end

    coeff = zeros(param.max_eigens, param.max_eigens, 4, param.num_bandits);

    stepno = 0;
    function step()
        a = find(cellfun(@length, p_hist) < 2, 1);
        p = task.get_pos();
        pt = p';
        if isempty(a)
            if rand() < param.epsilon
                a = randi([1, num_bandits], 1);
            else
                costs = zeros(1, param.num_bandits);
                [wx, wy] = ndgrid((1:param.max_eigens) * pi / param.chipwidth, (1:param.max_eigens) * pi / param.chipheight);
                for k = 1:param.num_bandits
                    pest = zeros(size(pt));
                    for n = 1:2:numel(p)
                        mdx = wx .* (-squeeze(coeff(:, :, 1, k)) .* sin(pt(n) .* wx) .* cos(pt(n + 1) .* wy) ...
                            +squeeze(coeff(:, :, 2, k)) .* cos(pt(n) .* wx) .* cos(pt(n + 1) .* wy) ...
                            -squeeze(coeff(:, :, 3, k)) .* sin(pt(n) .* wx) .* sin(pt(n + 1) .* wy) ...
                            +squeeze(coeff(:, :, 4, k)) .* cos(pt(n) .* wx) .* sin(pt(n + 1) .* wy));
                        dx = sum(sum(mdx));
                        mdy = wy .* (-squeeze(coeff(:, :, 1, k)) .* cos(pt(n) .* wx) .* sin(pt(n + 1) .* wy) ...
                            -squeeze(coeff(:, :, 2, k)) .* sin(pt(n) .* wx) .* sin(pt(n + 1) .* wy) ...
                            +squeeze(coeff(:, :, 3, k)) .* cos(pt(n) .* wx) .* cos(pt(n + 1) .* wy) ...
                            +squeeze(coeff(:, :, 4, k)) .* sin(pt(n) .* wx) .* cos(pt(n + 1) .* wy));
                        dy = sum(sum(mdy));
                        pest(n:n + 1) = pt(n:n + 1) + [dx, dy];
                    end
                    costs(k) = task.get_cost(pest');
                end
                [~, a] = min(costs);
            end
        end

        logging.log('sfb_pos', p);
        logging.log('sfb_actions', a);
        p_hist{a} = [p_hist{a}; reshape(pt, [], 1)];
        param.action(a);
        dp = task.get_pos() - p;
        logging.log('sfb_deltapos', dp);
        dp_hist{a} = [dp_hist{a}; reshape(dp', [], 1)];
        weights{a} = [weights{a} * param.sparse_decay; reshape(ones(size(dp')), [], 1)];
        l = length(dp_hist{a});
        if l >= 2
            if l > param.sparse_history
                p_hist{a} = p_hist{a}(end - param.sparse_history + 1:end, :);
                dp_hist{a} = dp_hist{a}(end - param.sparse_history + 1:end, :);
                weights{a} = weights{a}(end - param.sparse_history + 1:end, :);
            end
            b = dp_hist{a};
            pp = p_hist{a};
            A = zeros(length(b), param.max_eigens * param.max_eigens * 4);
            [wx, wy] = ndgrid((1:param.max_eigens) * pi / param.chipwidth, (1:param.max_eigens) * pi / param.chipheight);
            mdx = zeros(param.max_eigens, param.max_eigens, 4);
            mdy = zeros(param.max_eigens, param.max_eigens, 4);
            for n = 1:2:length(b)
                mdx(:, :, 1) = -wx .* sin(pp(n) .* wx) .* cos(pp(n + 1) .* wy);
                mdx(:, :, 2) = wx .* cos(pp(n) .* wx) .* cos(pp(n + 1) .* wy);
                mdx(:, :, 3) = -wx .* sin(pp(n) .* wx) .* sin(pp(n + 1) .* wy);
                mdx(:, :, 4) = wx .* cos(pp(n) .* wx) .* sin(pp(n + 1) .* wy);
                A(n, :) = reshape(mdx, 1, param.max_eigens * param.max_eigens * 4);
                mdy(:, :, 1) = -wy .* cos(pp(n) .* wx) .* sin(pp(n + 1) .* wy);
                mdy(:, :, 2) = -wy .* sin(pp(n) .* wx) .* sin(pp(n + 1) .* wy);
                mdy(:, :, 3) = wy .* cos(pp(n) .* wx) .* cos(pp(n + 1) .* wy);
                mdy(:, :, 4) = wy .* sin(pp(n) .* wx) .* cos(pp(n + 1) .* wy);
                A(n + 1, :) = reshape(mdy, 1, param.max_eigens * param.max_eigens * 4);
            end
            x = lasso(A, b, 'lambda', param.sparse_lambda, 'Weights', weights{a}); % TODO: add weights
            coeff(:, :, :, a) = reshape(x, param.max_eigens, param.max_eigens, 4, 1);
            if param.draw_sfb
                [wx, wy] = ndgrid((1:param.max_eigens) * pi, (1:param.max_eigens) * pi);
                p = zeros(size(imx));
                for i = 1:size(imx, 1)
                    for j = 1:size(imx, 2)
                        p(i, j) = sum(sum((squeeze(coeff(:, :, 1, a)) .* cos(imx(i, j) .* wx) .* cos(imy(i, j) .* wy) ...
                            +squeeze(coeff(:, :, 2, a)) .* sin(imx(i, j) .* wx) .* cos(imy(i, j) .* wy) ...
                            +squeeze(coeff(:, :, 3, a)) .* cos(imx(i, j) .* wx) .* sin(imy(i, j) .* wy) ...
                            +squeeze(coeff(:, :, 4, a)) .* sin(imx(i, j) .* wx) .* sin(imy(i, j) .* wy))));
                    end
                end
                minp = min(min(p));
                maxp = max(max(p));
                set(imobjs{a}, 'CData', (p - minp) * 255 / (maxp - minp));
                colormap hot;
            end
        end

        stepno = stepno + 1;
        if ~isempty(param.sparse_file) && param.sparse_save && mod(stepno, param.sparse_save_interval) == 0
            logging.message('Saving sfb model after %d steps', stepno);
            save(param.sparse_file, 'p_hist', 'dp_hist', 'weights');
        end
    end
end
