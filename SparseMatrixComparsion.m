n = 30;
myGuess = zeros(n * n, 1);
myTolerance = 0.00005;
myLambda = 1;
arrayIdx = 1;

ArcLengthStepSize = 0.1;
tic;
[myGuess0, myLambda0] = GuessInitialization(-0.1, 1, 1, 2 * pi * pi, 0.1, n);
U0 = fullNewtonFiniteElementMethod(myGuess0, myLambda0, myTolerance, n);
normUVec(arrayIdx) = norm(U0);
LambdaVec(arrayIdx) = myLambda0;
arrayIdx = arrayIdx + 1;

[myGuess1, myLambda1] = BetterGuessInitialization(U0, myLambda0, 0.1, n);
U1 = fullNewtonFiniteElementMethod(myGuess1, myLambda1, myTolerance, n);
normUVec(arrayIdx) = norm(U1);
LambdaVec(arrayIdx) = myLambda1;
arrayIdx = arrayIdx + 1;

[U2Guess, Lambda2Guess] = ALCGuess(U1, myLambda1, U0, myLambda0, ArcLengthStepSize, n);
S2 = sqrt((myLambda1 - myLambda0) * (myLambda1 - myLambda0) + norm(U1 - U0) * norm(U1 - U0)) + ArcLengthStepSize;
[U2Correct, Lambda2Correct] = fullNewtonAugmentedResidual(U2Guess, Lambda2Guess, S2, U0, myLambda0, 0, myTolerance, n);

Uk = U2Correct;
Lambdak = Lambda2Correct;

ArcLengthStepSize = 0.8;
while (Lambdak <= 60) && (Lambdak >= 0.5)
    [UkGuess, LambdakGuess] = ALCGuess(Uk, Lambdak, U0, myLambda0, ArcLengthStepSize, n);
    Sk = sqrt((Lambdak - myLambda0) * (Lambdak - myLambda0) + norm(Uk - U0) * norm(Uk - U0)) + ArcLengthStepSize;
    [Uk, Lambdak] = fullNewtonAugmentedResidual(UkGuess, LambdakGuess, Sk, U0, myLambda0, 0, myTolerance, n);
end

Time1 = toc;



ArcLengthStepSize = 0.1;
tic;
[myGuess0, myLambda0] = GuessInitialization(-0.1, 1, 1, 2 * pi * pi, 0.1, n);
U0 = fullNewtonSparseMatrix(myGuess0, myLambda0, myTolerance, n);
normUVec(arrayIdx) = norm(U0);
LambdaVec(arrayIdx) = myLambda0;
arrayIdx = arrayIdx + 1;

[myGuess1, myLambda1] = BetterGuessInitialization(U0, myLambda0, 0.1, n);
U1 = fullNewtonSparseMatrix(myGuess1, myLambda1, myTolerance, n);
normUVec(arrayIdx) = norm(U1);
LambdaVec(arrayIdx) = myLambda1;
arrayIdx = arrayIdx + 1;

[U2Guess, Lambda2Guess] = ALCGuess(U1, myLambda1, U0, myLambda0, ArcLengthStepSize, n);
S2 = sqrt((myLambda1 - myLambda0) * (myLambda1 - myLambda0) + norm(U1 - U0) * norm(U1 - U0)) + ArcLengthStepSize;
[U2Correct, Lambda2Correct] = fullNewtonAugmentedResidualSparseMatrix(U2Guess, Lambda2Guess, S2, U0, myLambda0, 0, myTolerance, n);

Uk = U2Correct;
Lambdak = Lambda2Correct;

ArcLengthStepSize = 0.8;
while (Lambdak <= 60) && (Lambdak >= 0.5)
    [UkGuess, LambdakGuess] = ALCGuess(Uk, Lambdak, U0, myLambda0, ArcLengthStepSize, n);
    Sk = sqrt((Lambdak - myLambda0) * (Lambdak - myLambda0) + norm(Uk - U0) * norm(Uk - U0)) + ArcLengthStepSize;
    [Uk, Lambdak] = fullNewtonAugmentedResidualSparseMatrix(UkGuess, LambdakGuess, Sk, U0, myLambda0, 0, myTolerance, n);
end

Time2 = toc;



