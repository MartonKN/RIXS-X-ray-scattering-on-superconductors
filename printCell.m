function printCell(c)
    mx = [];
    for i=1:length(c)
        mx = [mx, c{i}];
    end
    mx = mx(1:2,:),
end