function[U2Guess, Lambda2Guess] = ALCGuessSparseMatrix(iU1, iLambda1, iU0, iLambda0, ArcLengthStepSize, gridSize)
     h = 1 / (gridSize - 1);
    h2Inv = 1 / (h * h);
    variableNumber = gridSize * gridSize;
    J = eye(variableNumber, variableNumber);
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
                    J(row, col) = -4 * h2Inv + iLambda1 * (1 + 2 * iU1(i));
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
    
    % Construct ROLDri = deltaR(u1) / delta(lambda1):
    ROLDri = zeros(variableNumber, 1);
    for row = 1:gridSize
        for col = 1:gridSize
            idx = (row - 1) * gridSize + col;
            if atBorder(idx, gridSize)
                ROLDri(idx) = 0;
            else
                ROLDri(idx) = iU1(idx) * (1 + iU1(idx));
            end
        end
    end
   
    % Construct EOLDri = delta(eta)|1 / delta(lambda)|1 = -2 * (lambda1 - lambda0)
    EOLDri = -2 * (iLambda1 - iLambda0);
    
    % Construct EOUDri = delta(eta)|1 / delta(u)|1:
    EOUDri = -2 * (iU1 - iU0);
    
    % Construct Jhat:
    Jhat1FirstRow = [J ROLDri];
    Jhat1SecondRow = [transpose(EOUDri) EOLDri];
    Jhat1 = [Jhat1FirstRow; Jhat1SecondRow];
    Jhat1Sparse = sparse(Jhat1);
     
    % Construct EOSDri = delta(eta)|1 / delta(S)|1 = 2 * (S1 - S0) = +- 2 * sqrt((lambda1 - lambda0) * (lambda1 - lambda0) + norm(u1 - u0)):
    EOSDri = 2 * sqrt((iLambda1 - iLambda0) * (iLambda1 - iLambda0) + norm(iU1 - iU0) * norm(iU1 - iU0));
    
    % Construct RSVector = - [delta(R)|1 / delta(S) ; delta(Eta)|1 / delta(S)] 
    %                    = - [0 ; delta(Eta)|1 / delta(S)|1]
    %                    = - [0 ; 2 * (s1 - s0)]
    RSVector = zeros(variableNumber + 1, 1);
    RSVector(variableNumber + 1) = - EOSDri;
    
    % Calculate DriVector = [delta(u)|1 / delta(s)|1 ; delta(lambda)|1 / delta(s)|1]:
    DriVector = Jhat1Sparse \ RSVector;
    
    % Get UOSDri = delta(u)|1 / delta(s)|1;
    %     LOSDri = delta(lambda)|1 / delta(s)|1;
    % From DriVector;
    UOSDri = DriVector(1:variableNumber);
    LOSDri = DriVector(variableNumber + 1);
    
    U2Guess = iU1 + ArcLengthStepSize * UOSDri;
    Lambda2Guess = iLambda1 + ArcLengthStepSize * LOSDri;
end