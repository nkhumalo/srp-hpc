# README on WFN1

## Note on Monte-Carlo

In this code we use that when a calculation gets closer to the answer the
likelihood of guessing a step that improves the answer will decrease. If 
the likelihood becomes too small the range of the random numbers will be
shrunk. This increases the likelihood of guessing an improvement. A question
in this case is: how many successive failed guesses warrant reducing the 
random number range? To answer this question the same calculation was run
multiple times with numbers of successive fails. For each calculation the final
answer, the number of iterations to reach them, and the time to solution was
recorded. All calculations were performed at 0 temperature.

| fails | energy 1  | energy 2  | energy 3  | energy 4  | energy 5  | it 1  | it 2  | it 3  | it 4  | it 5  | time 1 | time 2 | time 3 | time 4 | time 5 |
| ----- | --------- | --------- | --------- | --------- | --------- | ----- | ----- | ----- | ----- | ----- | ------ | ------ | ------ | ------ | ------ |
| 20    | -1.104241 | -1.104241 | -1.104241 | -1.104241 | -1.104241 |    63 |    76 |    87 |    81 |    75 |    5.6 |    5.8 |    5.3 |    8.7 |    5.6 |
| 40    | -1.104241 | -1.104241 | -1.104241 | -1.104241 | -1.104241 |    85 |    69 |    75 |    75 |    70 |   11.1 |   15.4 |   10.1 |   10.4 |   10.8 |
| 60    | -1.104241 | -1.104241 | -1.104241 | -1.104241 | -1.104241 |    67 |    85 |    85 |    73 |    82 |   23.0 |   25.4 |   20.4 |   22.7 |   31.0 |
| 80    | -1.104241 | -1.104241 | -1.104241 | -1.104241 | -1.104241 |    89 |    74 |    67 |    68 |    74 |   21.5 |   20.6 |   20.0 |   21.1 |   26.7 |
| 100   | -1.104241 | -1.104241 | -1.104241 | -1.104241 | -1.104241 |    76 |    69 |    81 |    73 |    71 |   25.2 |   24.6 |   27.3 |   31.8 |   25.3 |

The results show that all calculations achieve the same energy independent
of the number of consecutive fails. The number of accepted steps to 
convergence on average was 76.4, 74.8, 78.4, 74.4, and 74.0 for 
consecutive fails ranging from 20 to 100. In other words the number
of steps to convergence is essentially independent of the number
of consecutive fails. However, as each failed step has to be retried
the time to solution is dependent on the number of consecutive fails,
generating average time solutions of 6.2, 11.6, 24.5, 22.0, and 26.8
seconds.

Repeating this experiment at a temperature of 0.2 starting from a 0 
temperature guess the answers are a bit different.

| fails | energy 1  | energy 2  | energy 3  | energy 4  | energy 5  | it 1  | it 2  | it 3  | it 4  | it 5  | time 1 | time 2 | time 3 | time 4 | time 5 |
| ----- | --------- | --------- | --------- | --------- | --------- | ----- | ----- | ----- | ----- | ----- | ------ | ------ | ------ | ------ | ------ |
| 20    | -1.337153 | -1.337125 | -1.337123 | -1.337123 | -1.337067 |   211 |   209 |   262 |   206 |   216 |    6.8 |    7.4 |   12.0 |    6.6 |    7.3 |
| 40    | -1.337153 | -1.337125 | -1.337125 | -1.337124 | -1.337122 |   227 |   188 |   179 |   249 |   240 |    8.5 |    8.3 |    6.3 |    7.5 |   11.1 |
