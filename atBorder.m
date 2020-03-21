function[atBorderAns] = atBorder(idx, gridSize)
    variableNumber = gridSize * gridSize;
    % outside the domain, the return value should be zero.
    if idx > variableNumber || idx < 1
        atBorderAns = 1;
        return;
    end
    % at the border, the U value should be zero.
    if mod(idx, gridSize) == 0 || mod(idx - 1, gridSize) == 0 || idx <= gridSize || (idx + gridSize) > variableNumber
        atBorderAns = 1;
        return;
    end
    atBorderAns = 0;
end