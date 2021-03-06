# Алгоритм глобального поиска
Результаты вычислительных эксперименов для оценки эффективности последовательного и параллельного алгоритмов поиска глобального минимума

#### Способы параллелизации
                
1. Разделение отрезка, на котором производится поиск глобального минимума, на подотрезки, внутри которых производится запуск последовательного алгоритма глобального поиска;
2. Распараллеливание вычисление внутренних характеристик последовательного алгоритма глобального поиска;
3. Распараллеливание и модифицирование алгоритма глобального поиска, обеспечивая одновременное выполнение нескольких испытаний.
                
                    
#### Координата х глобального минимума
                    
|№  | Real x | Sequential | MPI  1st way | MPI 2nd way | MPI  3rd way | TBB 1st way | TBB 3rd way |
|------------- | -------------|-------------|-------------|-------------|-------------|-------------|-------------|
|1  |10 | 9.99648 | 10.0012 | 9.99648 | 10.0018 | 10.0012 | 10.0018 |
|2  | 5.14575 | 5.14863 | 5.14671 | 5.14863 | 5.14557 | 5.14671 | 5.14557| 
|3  | -0.49139, -6.77458, 5.79179 | -6.78118 |  -0.48949 | -6.78118 | -0.49448 | -0.48949 | -0.49448|
|4  | 2.868  | 2.87484 |  2.8702 | 2.87484 | 2.86091 | 2.8702 | 2.86091|
|5  | 0.966086 | 0.9671 |  0.9671 | 0.9671 | 0.967406 | 0.9671 | 0.967406|
|6  | 0.679578 | 0.684807 |  0.688281 | 0.684807 | 0.684807 | 0.688281 | 0.684807|
|7  | 5.199776 | 5.19677 |  5.1967 | 5.19677 | 5.20088 | 5.1967 | 5.20088|
|8  | -0.80032, -7.08351, 5.48286 | 5.4838 |  5.48359 | 5.4838 | 5.48518 | 5.48359 | 5.48518|
|9  | 17.0392 | 17.0381|  17.0432 | 17.0381 | 17.0419 | 17.0432 | 17.0419|
|10  | 7.97867 | 7.98166 | 7.98245 | 7.98166 | 7.98143 | 7.98245 | 7.98143|
|11  | 2.09444, 4.18879  | 4.19226 |  4.18751 | 4.19226 | 4.18381 | 4.18751 | 4.18381|
|12  | 4.71239, 3.14159 | 4.71 |  3.14 | 4.71 | 4.71 | 3.14 | 4.71|
|13  | 0.70711 | 0.712251 | 0.712251 | 0.712251 | 0.71264 | 0.712251 | 0.71264|
|14  | 0.22488 | 0.219278 |  0.219278 | 0.219278 | 0.219278 | 0.219278 | 0.219278|
|15  | 2.41421 | 2.41259 |  2.42266 | 2.41259 | 2.41259 | 2.42266 | 2.41259|
|16 | 1.590721| 1.59762 |  1.59243 | 1.59762 | 1.59244 | 1.59243 | 1.59244|
|17  | -3, 3 | 3.0057 |  -2.99769 | 3.0057 | 2.9978 | -2.99769 | 2.9978|
|18  | 2 | 1.99946 |  2.00019 | 1.99946 | 2.00041 | 2.00019 | 2.00041|
|19  | 5.8728656 | 5.87662 |  5.87848 | 5.87662 | 5.87849 | 5.87848 | 5.87849|
|20  | 1.195137 | 1.19489 |  1.19277 | 1.19489 | 1.19903 | 1.19277 | 1.19903|

#### Время работы в секундах
                    
|№  | Sequential | MPI  1st way | MPI 2nd way | MPI  3rd way | TBB 1st way | TBB 3rd way |
|------------- |-------------|-------------|-------------|-------------|-------------|-------------|
|1  |1.98684 |  0.418129 | 0.0765019 | 0.746855  | 0.809128 | 9.52045 |
|2  | 0.351726 |  0.025808 | 0.0106301 | 0.0674728 | 0.0779781 | 0.856633 |
|3  | 1.98613 |  0.312087 | 0.101185  | 1.09666 | 0.969963 | 13.3693 |
|4  | 0.0580657 |  0.0109359| 0.0016188 | 0.0098205 | 0.0338061 | 0.205782 |
|5  | 0.0314042 |  0.0025829| 0.0011057 | 0.0057835| 0.0256853 | 0.118677 |
|6  | 2.12284 |  1.04983  | 0.0918233 | 0.762338 | 1.89965 | 11.9006 |
|7  | 0.515828 |  0.0325687 | 0.0138859 | 0.0785869 | 0.0914991 | 1.12353 |
|8  | 2.35078 |  0.28207 | 0.0889648  | 1.11156 | 0.994286 | 15.5373 |
|9  | 2.1886 |  0.31177 | 0.0824311 | 0.959739 | 0.945486 | 13.7701 |
|10  | 1.53283 |  0.184499 | 0.0317334 | 0.378982 | 0.488411   | 5.24417 |
|11  | 0.582684 |  0.0387124 | 0.0199062 | 0.133934 | 0.246798 | 1.9169 |
|12  | 0.434674 |  0.0321488  | 0.0133553  | 0.107945 | 0.155449 | 2.00538 |
|13  | 0.0175562 |  0.0020435| 0.000674 | 0.0026628| 0.0202261  | 0.0521011 |
|14  | 0.27692 |  0.01448  | 0.007001 | 0.0582406 | 0.0416089 | 0.926343 |
|15  | 1.06724 |  0.171527 | 0.0230642 | 0.266156 | 0.352577 | 4.1833|
|16 | 0.507303 |  0.0561889 | 0.0124977 | 0.103419 | 0.251011 | 1.45336 |
|17  | 0.69582 |  0.0723363 | 0.0200514  | 0.147168 | 0.341592 | 3.46478 |
|18  | 0.437395 |  0.0753891 | 0.0103295 | 0.0895955 | 0.214612 | 1.54216 |
|19  | 1.00528 |  0.0669877 | 0.019015  | 0.186115 | 0.163427 | 2.99287 |
|20  | 2.06142 |  0.987243 | 0.103779 | 0.875112 | 2.34018 | 11.3986 |
