n = 30;
sizeTimeIdx = 1;

for n = 15:35
    myGuess = zeros(n * n, 1);
    myTolerance = 0.00005;
    myLambda = 1;
    
    ArcLengthStepSize = 0.1;
    tic;
    [myGuess0, myLambda0] = GuessInitialization(-0.1, 1, 1, 2 * pi * pi, 0.1, n);
    U0 = fullNewtonSparseMatrix(myGuess0, myLambda0, myTolerance, n);
    
    [myGuess1, myLambda1] = BetterGuessInitializationSparseMatrix(U0, myLambda0, 0.1, n);
    U1 = fullNewtonSparseMatrix(myGuess1, myLambda1, myTolerance, n);
    
    [U2Guess, Lambda2Guess] = ALCGuessSparseMatrix(U1, myLambda1, U0, myLambda0, ArcLengthStepSize, n);
    S2 = sqrt((myLambda1 - myLambda0) * (myLambda1 - myLambda0) + norm(U1 - U0) * norm(U1 - U0)) + ArcLengthStepSize;
    [U2Correct, Lambda2Correct] = fullNewtonAugmentedResidualSparseMatrix(U2Guess, Lambda2Guess, S2, U0, myLambda0, 0, myTolerance, n);
    
    Uk = U2Correct;
    Lambdak = Lambda2Correct;
    
    ArcLengthStepSize = 0.8;
    while (Lambdak <= 60) && (Lambdak >= 0.5)
        [UkGuess, LambdakGuess] = ALCGuessSparseMatrix(Uk, Lambdak, U0, myLambda0, ArcLengthStepSize, n);
        Sk = sqrt((Lambdak - myLambda0) * (Lambdak - myLambda0) + norm(Uk - U0) * norm(Uk - U0)) + ArcLengthStepSize;
        [Uk, Lambdak] = fullNewtonAugmentedResidualSparseMatrix(UkGuess, LambdakGuess, Sk, U0, myLambda0, 0, myTolerance, n);
    end

    timeVec(sizeTimeIdx) = toc;
    sizeVec(sizeTimeIdx) = n;
    sizeTimeIdx = sizeTimeIdx + 1;
end

sizeVec = transpose(sizeVec);
timeVec = transpose(timeVec);

[xDim, yDim] = size(sizeVec);
for idx = 1:xDim
    sizeVec(idx) = sizeVec(idx) * sizeVec(idx) * sizeVec(idx);
end

plot(sizeVec, timeVec, '.');

myTable = table(sizeVec, timeVec);
filename = 'sparseMatrixTimeSize.xlsx';
writetable(myTable, filename, 'Sheet', 1, 'Range', 'A1');