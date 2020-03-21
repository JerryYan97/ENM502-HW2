function[GuessU1, GuessLambda1] = BetterGuessInitialization(iGuessU0, iGuessLambda0, iLambdaStepSize, gridSize) % AC process.
    GuessLambda1 = iGuessLambda0 + iLambdaStepSize;
    % Construct J0:
    variableNumber = gridSize * gridSize;
    J0 = eye(variableNumber, variableNumber);
    h = 1 / (gridSize - 1);
    h2Inv = 1 / (h * h);
    for row = 1:variableNumber
        for col = 1:variableNumber
            i = row;
            j = col;
            if atBorder(i, gridSize)
                if i == j
                    J0(row, col) = 1;
                else
                    J0(row, col) = 0;
                end
            else
                if j == (i - gridSize)
                    % case a
                    J0(row, col) = h2Inv;
                elseif j == (i - 1)
                    % case b
                    J0(row, col) = h2Inv;
                elseif j == i
                    % case c
                    J0(row, col) = -4 * h2Inv + iGuessLambda0 * (1 + 2 * iGuessU0(i));
                elseif j == i + 1
                    % case d
                    J0(row, col) = h2Inv;
                elseif j == i + gridSize
                    % case e
                    J0(row, col) = h2Inv;
                else
                    J0(row, col) = 0;
                end
            end
        end
    end
    
    % Construct ROLDri = - deltaR(u0) / delta(lambda0):
    ROLDri = zeros(variableNumber, 1);
    for row = 1:gridSize
        for col = 1:gridSize
            idx = (row - 1) * gridSize + col;
            if atBorder(idx, gridSize)
                ROLDri(idx) = 0;
            else
                ROLDri(idx) = - iGuessU0(idx) * (1 + iGuessU0(idx));
            end
        end
    end
    
    % Get UOLDri = delta(u0) / delta(lambda0):
    UOLDri = J0 \ ROLDri;
    
    % Get U1 first guess:
    GuessU1 = iGuessU0 + UOLDri * iLambdaStepSize;
    
    % Get U1 convergence:
    GuessU1 = fullNewtonFiniteElementMethod(GuessU1, GuessLambda1, 0.00005, gridSize);
end