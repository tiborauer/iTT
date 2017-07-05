function out = check_licence(func)
for i = 1:numel(func)
    try
        eval(func{i});
    catch em
        switch em.identifier
            case 'MATLAB:license:checkouterror'
                out(i) = 1;
                return
        end
    end
    out(i) = 0;
end