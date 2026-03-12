# Class 5- Data Viz with ggplot
Joshua Khalil (PID:A17784122)

- [Background](#background)
- [Gene Expression Plot](#gene-expression-plot)
- [Going Further](#going-further)
- [Custom Plot](#custom-plot)

## Background

There are lots of ways to make figures in R. These include so called
“base R” graphics (e.g. `plot()`) qnd tons of add on packages like
**ggplot2**.

For example here, we make the same plot with both:

``` r
head(cars)
```

      speed dist
    1     4    2
    2     4   10
    3     7    4
    4     7   22
    5     8   16
    6     9   10

``` r
plot(cars)
```

![](Class05_files/figure-commonmark/unnamed-chunk-2-1.png)

First I need to install the package with the command
`install.packages()`

> **N.B** We never run an install cmd in quarto code chunk or we will
> end up re-installing packages many times which is not what we want!

Everytime we want to use one of these “add-on” packages we need to load
it up in R with the `library()` function:

``` r
library(ggplot2)
```

``` r
ggplot(cars)
```

![](Class05_files/figure-commonmark/unnamed-chunk-4-1.png)

Every ggplot needs at least three things:

- The **data**, the information you want plotted
- The **aes**thitics, how the data map to the plot
- The **geom**etry, the type of plot

``` r
p <- ggplot(cars) +
  aes(x=speed, y=dist) +
 geom_point() +
geom_smooth(method="lm", se= FALSE) +
  labs(title="Stopping distance of old cars",
       subtitle = "Data from the `cars` object",
       x= "Speed (mph)",
       y= "Distance (ft)")
```

Add a line to better show relationship

Render it out

``` r
p + theme_bw()
```

    `geom_smooth()` using formula = 'y ~ x'

![](Class05_files/figure-commonmark/unnamed-chunk-6-1.png)

## Gene Expression Plot

We can read the input data from ther class website

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

A first version plot

``` r
ggplot(genes) +
  aes(Condition1, Condition2 ) +
  geom_point()
```

![](Class05_files/figure-commonmark/unnamed-chunk-8-1.png)

``` r
table( genes$State)
```


          down unchanging         up 
            72       4997        127 

Version 2 lets color by `State` so we can see the “up” and “down”
significant genes compared to all the unchanging genes

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col=State) +
  geom_point()
```

![](Class05_files/figure-commonmark/unnamed-chunk-10-1.png)

Version 3 plot, let’s modify the default colors to our preferences

``` r
ggplot(genes) +
  aes(Condition1, Condition2, col=State) +
  geom_point() + 
  scale_colour_manual( values=c("Red","Grey","yellow") ) +
labs(x="Control (no drug)",
      y="Drug",
    title="Gene Expression Changes upon GLP-1 drug" ) +
  theme_bw()
```

![](Class05_files/figure-commonmark/unnamed-chunk-11-1.png)

## Going Further

Let’s have a look at the famous **gapminder** data set

``` r
# File location online
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
```

``` r
head(gapminder,3)
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007

Version 1

``` r
ggplot(gapminder)+
  aes(gdpPercap, lifeExp, col=continent) +
  geom_point(alpha=0.3,)
```

![](Class05_files/figure-commonmark/unnamed-chunk-14-1.png)

Let’s “facet” (i.e. make a separate plot) by continent rather than the
big hot mess

``` r
ggplot(gapminder)+
  aes(gdpPercap, lifeExp, col=continent) +
  geom_point(alpha=0.3,) +
  facet_wrap(~continent)
```

![](Class05_files/figure-commonmark/unnamed-chunk-15-1.png)

## Custom Plot

How big is this gapminder dataset?

``` r
nrow(gapminder)
```

    [1] 1704

``` r
ncol(gapminder)
```

    [1] 6

I want to “filter” down to a subset of this data. I will use the
**dplyr** package to help me.

First, I need to install it and load it up… `install.packages("dplyr")`
and then `library(dplyr)`

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
gapminder_2007 <- filter(gapminder, year ==2007)
head(gapminder_2007)
```

          country continent year lifeExp      pop  gdpPercap
    1 Afghanistan      Asia 2007  43.828 31889923   974.5803
    2     Albania    Europe 2007  76.423  3600523  5937.0295
    3     Algeria    Africa 2007  72.301 33333216  6223.3675
    4      Angola    Africa 2007  42.731 12420476  4797.2313
    5   Argentina  Americas 2007  75.320 40301927 12779.3796
    6   Australia   Oceania 2007  81.235 20434176 34435.3674

``` r
filter(gapminder_2007, country=="Ireland")
```

      country continent year lifeExp     pop gdpPercap
    1 Ireland    Europe 2007  78.885 4109086     40676

``` r
filter(gapminder, year==2007, country=="Ireland")
```

      country continent year lifeExp     pop gdpPercap
    1 Ireland    Europe 2007  78.885 4109086     40676

``` r
filter(gapminder, year==1977, country=="United States")
```

            country continent year lifeExp       pop gdpPercap
    1 United States  Americas 1977   73.38 220239000  24072.63

> Q. make a plot comparing 1977 and 2007 for all countries

``` r
plotss <- filter(gapminder, year %in% c(1977,2007))
ggplot(plotss)+
  aes(gdpPercap, lifeExp, col=continent) +
  geom_point() +
  facet_wrap(~year)
```

![](Class05_files/figure-commonmark/unnamed-chunk-21-1.png)
