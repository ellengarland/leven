# leven

Levenshtein distance calculations

## Installation
You should be able to install this by installing the devtools package
```
install.packages("devtools")
```

then compile and install the leven package:
```
library("devtools")
install_github("ellengarland/leven")
```

## Usage
Once installed, add it to your workspace with
```
library("leven")
```
and then use it:

```
> leven(c("k", "i", "t", "t", "e", "n"), c("s", "i", "t", "t", "i", "n", "g"))
[1] 3
```
