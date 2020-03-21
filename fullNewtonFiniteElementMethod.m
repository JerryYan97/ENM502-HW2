function[fullNewTonAns] = fullNewtonFiniteElementMethod(initGuess, lambda,tolerence, gridSize)
    h = 1 / (gridSize - 1);
    h2Inv = 1 / (h * h);
    variableNumber = gridSize * gridSize;
    J = eye(variableNumber, variableNumber);
    error = 10000;
    deltaU = zeros(variableNumber, 1);
    R = zeros(variableNumber, 1);
    U = initGuess;
    % init U: make the U value at border become zero.
    for idx = 1:variableNumber
        U(idx) = uValue(U, idx, gridSize);
    end
    
    while error > tolerence
        % setup J matrix:
        % the idx here is for a element in the matrix
        for row = 1:variableNumber
            for col = 1:variableNumber
                idx = (row - 1) * gridSize + col;
                i = row;
                j = col;
                if atBorder(i, gridSize)
                    if i == j
                        J(row, col) = 1;
                    else
                        J(row, col) = 0;
                    end
                else
                    if j == (i - gridSize)
                        % case a
                        J(row, col) = h2Inv;
                    elseif j == (i - 1)
                        % case b
                        J(row, col) = h2Inv;
                    elseif j == i
                        % case c
                        J(row, col) = -4 * h2Inv + lambda * (1 + 2 * U(i));
                    elseif j == i + 1
                        % case d
                        J(row, col) = h2Inv;
                    elseif j == i + gridSize
                        % case e
                        J(row, col) = h2Inv;
                    else
                        J(row, col) = 0;
                    end
                end
            end
        end
        % setup -R(uk) vector:
        % the idx here is for a point in the grid
        for row = 1:gridSize
            for col = 1:gridSize
                idx = (row - 1) * gridSize + col;
                if atBorder(idx, gridSize)
                    R(idx) = -U(idx);
                else
                    idxMin1 = idx - 1;
                    idxAdd1 = idx + 1;
                    idxMinNx = idx - gridSize;
                    idxAddNx = idx + gridSize;
                    UiAdd1 = U(idxAdd1);
                    Ui = U(idx);
                    UiMin1 = U(idxMin1);
                    UiAddNx = U(idxAddNx);
                    UiMinNx = U(idxMinNx);
                    R(idx) = -((UiAdd1 - 2 * Ui + UiMin1) * h2Inv + (UiAddNx - 2 * Ui + UiMinNx) * h2Inv + lambda * Ui * (1 + Ui));
                end                
            end
        end
        % solve deltaU:
        deltaU = J \ R; 
        
        % update U:
        U = U + deltaU;
        % make the U value at border become zero.
        for idx = 1:variableNumber
            U(idx) = uValue(U, idx, gridSize);
        end
        
        % update error:
        errVec = zeros((gridSize - 2) * (gridSize - 2), 1);
        errVecIdx = 1;
        for row = 2:(gridSize - 1)
            for col = 2:(gridSize - 1)
                idx = (row - 1) * gridSize + col;
                errVec(errVecIdx) = deltaU(idx);
                errVecIdx = errVecIdx + 1;
            end
        end
        error = norm(errVec);
    end
    for idx = 1:variableNumber
        U(idx) = uValue(U, idx, gridSize);
    end
    fullNewTonAns = U;
end