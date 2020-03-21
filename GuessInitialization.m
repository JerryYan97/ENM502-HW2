function[GuessU0, GuessLambda0] = GuessInitialization(A, m, n, LambdaSol, lambdaOffset,gridSize)
    % GuessAlmbda0 = 2 * pi * pi + lambdaOffset;
    GuessLambda0 = LambdaSol + lambdaOffset;
    GuessU0 = zeros(gridSize * gridSize, 1);
    for row = 1:gridSize
        for col = 1:gridSize
            idx = (row - 1) * gridSize + col;
            x = (col - 1) / (gridSize - 1);
            y = (row - 1) / (gridSize - 1);
            GuessU0(idx) = A * sin(m * pi * x) * sin(n * pi * y);
        end
    end
end