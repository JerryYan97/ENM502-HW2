gridSize = 30;
myGuess = zeros(gridSize * gridSize, 1);
myTolerance = 0.00005;
myLambda = 1;
ArcLengthStepSize = 0.1;
m = 1;
n = 1;
arrayIdx = 1;

[myGuess0, myLambda0] = GuessInitialization(0.1, m, n, 2 * pi * pi, -0.1, gridSize);
U0 = fullNewtonFiniteElementMethod(myGuess0, myLambda0, myTolerance, gridSize);
[myGuess1, myLambda1] = BetterGuessInitialization(U0, myLambda0, -0.1, gridSize);
U1 = fullNewtonFiniteElementMethod(myGuess1, myLambda1, myTolerance, gridSize);
[U2Guess, Lambda2Guess] = ALCGuess(U1, myLambda1, U0, myLambda0, ArcLengthStepSize, gridSize);
S2 = sqrt((myLambda1 - myLambda0) * (myLambda1 - myLambda0) + norm(U1 - U0) * norm(U1 - U0)) + ArcLengthStepSize;
[U2Correct, Lambda2Correct] = fullNewtonAugmentedResidual(U2Guess, Lambda2Guess, S2, U0, myLambda0, 0, myTolerance, gridSize);

Uk = U2Correct;
Lambdak = Lambda2Correct;

ArcLengthStepSize = 0.8;
while (Lambdak <= 60) && (Lambdak >= 1)
    [UkGuess, LambdakGuess] = ALCGuess(Uk, Lambdak, U0, myLambda0, ArcLengthStepSize, gridSize);
    Sk = sqrt((Lambdak - myLambda0) * (Lambdak - myLambda0) + norm(Uk - U0) * norm(Uk - U0)) + ArcLengthStepSize;
    [Uk, Lambdak] = fullNewtonAugmentedResidual(UkGuess, LambdakGuess, Sk, U0, myLambda0, 0, myTolerance, gridSize);
    arrayIdx = arrayIdx + 1;
end

U = zeros(gridSize, gridSize);
for row = 1:gridSize
    for col = 1:gridSize
        idx = (row - 1) * gridSize + col;
        U(row, col) = Uk(idx);
    end
end

figure(5)
contourf(U)
colorbar
title(['Positive. m = ' num2str(m) ', n = ' num2str(n) '. Lambda = ' num2str(Lambdak) '.'])