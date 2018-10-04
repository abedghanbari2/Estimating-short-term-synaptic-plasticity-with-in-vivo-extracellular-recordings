function y = logitfxn(x)
    y = log(x ./ (1-x));