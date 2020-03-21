function[uReturn] = uValue(U, idx, gridSize)
% return U value if it is in the domain.
% at the border, the U value should be zero.
% outside the domain, the return value should be zero.
    variableNumber = gridSize * gridSize;
    % outside the domain, the return value should be zero.
    if idx > variableNumber || idx < 1
        uReturn = 0;
        return;
    end
    % at the border, the U value should be zero.
    if mod(idx, gridSize) == 0 || mod(idx - 1, gridSize) == 0 || idx <= gridSize || (idx + gridSize) > variableNumber
        uReturn = 0;
        return;
    end
    % return U value if it is in the domain.
    uReturn = U(idx);
end