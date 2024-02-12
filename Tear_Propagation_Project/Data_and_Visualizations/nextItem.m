function line = nextItem(fileID, item)
    n = length(item);
    switch nargout
        case 1
            tline = fgetl(fileID);
            %lineCounter = 1;
            while ischar(tline)
                %disp(tline)
                if strncmp(tline, item, n)
                    line = tline;
                    break;
                end
                % Read next line
                tline = fgetl(fileID);
                %lineCounter = lineCounter + 1;
            end
        otherwise
            tline = fgetl(fileID);
            %lineCounter = 1;
            while ischar(tline)
                %disp(tline)
                if strncmp(tline, item, n)
                    break;
                end
                % Read next line
                tline = fgetl(fileID);
                %lineCounter = lineCounter + 1;
            end
            %ftell(fileID)
    end
end