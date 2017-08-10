function outputbw = xihuaquzhi(inputbw, Drection)

outputbw = inputbw;
[hbw,wbw] = size(inputbw);

if(strcmp(Drection,'vd&up'))
    for j = 1:wbw
        if sum(inputbw(:,j))>1
            for i = 1:hbw
                if inputbw(i,j) == 1
                    inputbw(:,j) = 0;
                    inputbw(i,j) = 1;
                    outputbw(:,j) = inputbw(:,j);
                end
            end
        end
    end
end
if(strcmp(Drection,'vd&down'))
    for j = 1:wbw
        if sum(inputbw(:,j))>1
            for i = hbw:-1:1
                if inputbw(i,j) == 1
                    inputbw(:,j) = 0;
                    inputbw(i,j) = 1;
                    outputbw(:,j) = inputbw(:,j);
                end
            end
        end
    end
end
if(strcmp(Drection,'hd&left'))
    for i = 1:hbw
        if sum(inputbw(i,:))>1
            for j = 1:wbw
                if inputbw(i,j) == 1
                    inputbw(i,:) = 0;
                    inputbw(i,j) = 1;
                    outputbw(i,:) = inputbw(i,:);
                end
            end
        end
    end
end
if(strcmp(Drection,'hd&right'))
    for i = 1:hbw
        if sum(inputbw(i,:))>1
            for j = wbw:-1:1
                if inputbw(i,j) == 1
                    inputbw(i,:) = 0;
                    inputbw(i,j) = 1;
                    outputbw(i,:) = inputbw(i,:);
                end
            end
        end
    end
end
end

