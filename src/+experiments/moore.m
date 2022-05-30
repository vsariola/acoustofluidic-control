function moore(varargin)
    order = 5;
    width = 7; % mm
    height = 6;
    margin = 0.5; % mm from the edges

    a = 1 + 1i;
    b = 1 - 1i;

    % Generate point sequence, based on
    % https://blogs.mathworks.com/steve/2012/01/25/generating-hilbert-curves/
    z = 0;
    for k = 1:order - 1
        w = 1i * conj(z);
        z = [w - a; z - b; z + a; b - w] / 2;
    end

    % Moore curve: four Hilbert curves in the quadrants, returning to
    % original point
    z = [z - 1 + 1i; z + 1 + 1i; -z + 1 - 1i; -z - 1 - 1i] / 2;

    paths = ([real(z) imag(z)] + 1) / 2 .* ([width height] - 2 * margin) + margin;

    control_loop('logname', 'moore', ...
        'paths', paths, ...
        'draw', true, ...
        'controller', 'sfb', ...
        'sparse_file', 'sfb_model.mat', ...
        'sparse_load', false, ...
        'sparse_save', true, ...
        'save_images', true, ...
        varargin{:});
end
