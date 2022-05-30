function ret = flat_model(varargin)
    import logging.*;

    parser = inputParser;
    parser.addParameter('amp', 0.032);
    parser.addParameter('minfreq', 65e3);
    parser.addParameter('maxfreq', 700e3);
    parser.addParameter('numfreq', 100);
    parser.addParameter('chipwidth', 7);
    parser.addParameter('chipheight', 6);
    parser.KeepUnmatched = true;
    parse(parser, varargin{:});
    r = parser.Results;
    logging.message('%s\n%s', mfilename, third_party.struct2str(r));

    amps = ones(r.numfreq, 1) * r.amp;
    freqs = linspace(r.minfreq, r.maxfreq, r.numfreq);

    ret = struct( ...
        'chipwidth', r.chipwidth, ...
        'chipheight', r.chipheight, ...
        'numfreq', r.numfreq, ...
        'freq', freqs, ...
        'amp', amps);
end
