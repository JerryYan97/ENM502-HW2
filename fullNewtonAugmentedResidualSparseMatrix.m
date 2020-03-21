function[U2Correct, Lambda2Correct] = fullNewtonAugmentedResidualSparseMatrix(iU2Guess, iLambda2Guess, S2, iU0, iLambda0, S0, tolerance, gridSize)
    h = 1 / (gridSize - 1);
    h2Inv = 1 / (h * h);
    variableNumber = gridSize * gridSize;
    J = eye(variableNumber, variableNumber);
    error = 10000;
    U = iU2Guess;
    Lambda = iLambda2Guess;
    while error > tolerance
        % init U: make the U value at border become zero.
        for idx = 1:variableNumber
            U(idx) = uValue(U, idx, gridSize);
        end
        
        % setup J matrix:
        % the idx here is for a element in the matrix
        for row = 1:variableNumber
            for col = 1:variableNumber
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
                        J(row, col) = -4 * h2Inv + Lambda * (1 + 2 * U(i));
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
        
        % Construct ROLDri = deltaR(u)|k / delta(lambda)|k:
        ROLDri = zeros(variableNumber, 1);
        for row = 1:gridSize
            for col = 1:gridSize
                idx = (row - 1) * gridSize + col;
                if atBorder(idx, gridSize)
                    ROLDri(idx) = 0;
                else
                    ROLDri(idx) = U(idx) * (1 + U(idx));
                end
            end
        end
        
        % Construct EOLDri = delta(eta)|k / delta(lambda)|k = -2 * (lambdak - lambda0)
        EOLDri = -2 * (Lambda - iLambda0);
        
        % Construct EOUDri = delta(eta)|k / delta(u)|0:
        EOUDri = -2 * (U - iU0);
        
        % Construct Jhat:
        JhatFirstRow = [J ROLDri];
        JhatSecondRow = [transpose(EOUDri) EOLDri];
        Jhat = [JhatFirstRow; JhatSecondRow];
        JhatSparse = sparse(Jhat);
        
        % Construct REtaKVector = - [R|k; Eta|k];
        % Construct RkVector = R|k;
        RkVector = zeros(variableNumber, 1);
        for row = 1:gridSize
            for col = 1:gridSize
                idx = (row - 1) * gridSize + col;
                if atBorder(idx, gridSize)
                    RkVector(idx) = U(idx);
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
                    RkVector(idx) = ((UiAdd1 - 2 * Ui + UiMin1) * h2Inv + (UiAddNx - 2 * Ui + UiMinNx) * h2Inv + Lambda * Ui * (1 + Ui));
                end
            end
        end
        
        % Construct Etak = Eta|k;
        Etak = (S2 - S0) * (S2 - S0) - (Lambda - iLambda0) * (Lambda - iLambda0) - norm(U - iU0) * norm(U - iU0);
        
        REtaKVector = - [RkVector; Etak];
        
        % Calculate delta(U)|k and delta(Lambda)|k:
        % DeltaULambdaVec = [delta(U)|k; delta(Lambda)|k];
        DeltaULambdaVec = JhatSparse \ REtaKVector;
        DeltaU = DeltaULambdaVec(1:variableNumber);
        DeltaLambda = DeltaULambdaVec(variableNumber + 1);
        
        U = U + DeltaU;
        Lambda = Lambda + DeltaLambda;
        
        error = norm(DeltaULambdaVec);
    end
    for idx = 1:variableNumber
        U(idx) = uValue(U, idx, gridSize);
    end
    U2Correct = U;
    Lambda2Correct = Lambda;
end