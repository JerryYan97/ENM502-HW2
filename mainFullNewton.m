n = 30;
myGuess = zeros(n * n, 1);
myTolerance = 0.00005;
myLambda = 1;
ArcLengthStepSize = 0.1;
[myGuess0, myLambda0] = GuessInitialization(-0.1, 1, 2, 5 * pi * pi, 0.1, n);
U0 = fullNewtonFiniteElementMethod(myGuess0, myLambda0, myTolerance, n);
[myGuess1, myLambda1] = BetterGuessInitialization(myGuess0, myLambda0, 0.1, n);
U1 = fullNewtonFiniteElementMethod(myGuess1, myLambda1, myTolerance, n);
[U2Guess, Lambda2Guess] = ALCGuess(U1, myLambda1, U0, myLambda0, ArcLengthStepSize, n);
S2 = sqrt((myLambda1 - myLambda0) * (myLambda1 - myLambda0) + norm(U1 - U0) * norm(U1 - U0)) + ArcLengthStepSize;
[U2Correct, Lambda2Correct] = fullNewtonAugmentedResidual(U2Guess, Lambda2Guess, S2, U0, myLambda0, 0, myTolerance, n);

U = zeros(n, n);
for row = 1:n
    for col = 1:n
        idx = (row - 1) * n + col;
        U(row, col) = U2Correct(idx);
    end
end

figure(5)
contourf(U)
colorbar
title(['Contour map for 30 by 30 grid, lambda is 14.9404'])
