function [ recon_params ] = yarra_read_mode_section(mode_path, section)
%% Reads a single section from an ini file, doing simple data conversions
%  Strings can be converted to numbers successfully with str2num are
%  converted
section_values = inifile(mode_path,'readall');
recon_params = struct;
for i=1:size(section_values)-1
    if ((section_values{i,1}) == string(lower(section)))
        [num,status] = str2num(section_values{i,4}); %#ok<ST2NM>
        if (status == 1)
            recon_params.(section_values{i,3}) = num;
        else
            recon_params.(section_values{i,3}) = section_values{i,4};
        end
    end
end
end