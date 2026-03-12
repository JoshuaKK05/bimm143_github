# Class 18: Pertussis Mini Project
Joshua Khalil (PID A17784122)

## Background

Pertussis whooping cough is a common lung infection caused by the
bacteria B. Pertussis. It can infect anyone but it is most deadly for
infants (under 1 year of age)

## CDC Tracking data

The CDC tracks the number of Pertussis cases:

> Q1. With the help of the R “addin” package datapasta assign the CDC
> pertussis case number data to a data frame called cdc and use ggplot
> to make a plot of cases numbers over time.

I want a plot of year vs cases

``` r
library(ggplot2)
ggplot(cdc) +
  aes(year, cases) +
  geom_point() +
geom_line()
```

![](Class18_files/figure-commonmark/unnamed-chunk-2-1.png)

> Q2. Using the ggplot geom_vline() function add lines to your previous
> plot for the 1946 introduction of the wP vaccine and the 1996 switch
> to aP vaccine (see example in the hint below). What do you notice?
> (add annotation lines for the major milestones of wP vaccination
> roll-out (1946) and the switch to the aP vaccine (1996).)

``` r
library(ggplot2)
ggplot(cdc) +
  aes(year, cases) +
  geom_point() +
geom_line() +
geom_vline(xintercept=1946, col= "blue", lty=2)  +
geom_vline(xintercept= 1996, col = "red", lty=2)+
geom_vline(xintercept=2020, col="grey", lty=2)
```

![](Class18_files/figure-commonmark/unnamed-chunk-3-1.png)

I notice that cases are decreasing from the 1946 point from over 100,000
cases, down to below 10,000 cases, after the development of the wP
vaccination that significantly decreased cases from before 1946 to after
1946… However as time progresses the number of cases begins to rise.

> Q3. Describe what happened after the introduction of the aP vaccine?
> Do you have a possible explanation for the observed trend?

After the introduction of the aP vaccines, the cases began to rise again
up to 20,000 cases. Once 2020 hit, the cases went back down to less than
10,000 cases. This could be because of resistance of the bacteria
adapting and mutating over time, or it could be because of less people
taking vaccinations out of fear. Or the aP vaccine could not be as
strong as the old one, the aP vaccine based on the evidence suggests
that the effects of the aP vaccination have much less of an effect on
pertussis than the original wP vaccination, proving to be less
effective.

## Exploring CMI-PB data

The CMI-PB Project \< https://www.cmi-pb.org/\> mission is to provide
the scientific community with a comprehensive, high-quality and freely
accessible resource of Pertussis booster vaccination.

Basically make available a large data set on the immune response to
Pertussis. They use a booster vaccination as a proxy for Pertussis
vaccination.

They make their data available as JSON format API. We can read this into
R with the `read_json()` function from **jsonlite** package

``` r
library(jsonlite)
subject <- read_json("https://www.cmi-pb.org/api/v5_1/subject", simplifyVector = TRUE)
head(subject)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          2          wP         Female Not Hispanic or Latino White
    3          3          wP         Female                Unknown White
    4          4          wP           Male Not Hispanic or Latino Asian
    5          5          wP           Male Not Hispanic or Latino Asian
    6          6          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset
    1    1986-01-01    2016-09-12 2020_dataset
    2    1968-01-01    2019-01-28 2020_dataset
    3    1983-01-01    2016-10-10 2020_dataset
    4    1988-01-01    2016-08-29 2020_dataset
    5    1991-01-01    2016-08-29 2020_dataset
    6    1988-01-01    2016-10-10 2020_dataset

> Q4. How many aP and wP infancy vaccinated subjects are in the dataset?

``` r
table(subject$infancy_vac)
```


    aP wP 
    87 85 

> Q5. How many Male and Female subjects/patients are in the dataset?

``` r
table(subject$biological_sex)
```


    Female   Male 
       112     60 

> Q6. What is the breakdown of race and biological sex (e.g. number of
> Asian females, White males etc…)? Is this representative of the US
> Population

No it is not representative of the US population

``` r
table(subject$race, subject$biological_sex)
```

                                               
                                                Female Male
      American Indian/Alaska Native                  0    1
      Asian                                         32   12
      Black or African American                      2    3
      More Than One Race                            15    4
      Native Hawaiian or Other Pacific Islander      1    1
      Unknown or Not Reported                       14    7
      White                                         48   32

We can now read more tables from the CMI-PB database

``` r
specimen <- read_json("https://www.cmi-pb.org/api/v5_1/specimen",
                      simplifyVector = TRUE)
ab_titer <- read_json("https://www.cmi-pb.org/api/v5_1/plasma_ab_titer",
                      simplifyVector = TRUE)
```

``` r
head(specimen)
```

      specimen_id subject_id actual_day_relative_to_boost
    1           1          1                           -3
    2           2          1                            1
    3           3          1                            3
    4           4          1                            7
    5           5          1                           11
    6           6          1                           32
      planned_day_relative_to_boost specimen_type visit
    1                             0         Blood     1
    2                             1         Blood     2
    3                             3         Blood     3
    4                             7         Blood     4
    5                            14         Blood     5
    6                            30         Blood     6

``` r
head(ab_titer)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgE               FALSE   Total 1110.21154       2.493425
    2           1     IgE               FALSE   Total 2708.91616       2.493425
    3           1     IgG                TRUE      PT   68.56614       3.736992
    4           1     IgG                TRUE     PRN  332.12718       2.602350
    5           1     IgG                TRUE     FHA 1887.12263      34.050956
    6           1     IgE                TRUE     ACT    0.10000       1.000000
       unit lower_limit_of_detection
    1 UG/ML                 2.096133
    2 IU/ML                29.170000
    3 IU/ML                 0.530000
    4 IU/ML                 6.205949
    5 IU/ML                 4.679535
    6 IU/ML                 2.816431

To make sense of all this data we need to “join” the tables together.
Only then will you know that a given antibody measurement (from the
`ab_titer` table) was collected on certain days (from the `specimen`
table) from a certain wP or aP subject (from the `subject` table).

We can use **dplyr** and the `*join()` family function to do this.

> Q9. Complete the code to join specimen and subject tables to make a
> new merged data frame containing all specimen records along with their
> associated subject details:

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
meta<- inner_join(subject,specimen)
```

    Joining with `by = join_by(subject_id)`

``` r
head(meta)
```

      subject_id infancy_vac biological_sex              ethnicity  race
    1          1          wP         Female Not Hispanic or Latino White
    2          1          wP         Female Not Hispanic or Latino White
    3          1          wP         Female Not Hispanic or Latino White
    4          1          wP         Female Not Hispanic or Latino White
    5          1          wP         Female Not Hispanic or Latino White
    6          1          wP         Female Not Hispanic or Latino White
      year_of_birth date_of_boost      dataset specimen_id
    1    1986-01-01    2016-09-12 2020_dataset           1
    2    1986-01-01    2016-09-12 2020_dataset           2
    3    1986-01-01    2016-09-12 2020_dataset           3
    4    1986-01-01    2016-09-12 2020_dataset           4
    5    1986-01-01    2016-09-12 2020_dataset           5
    6    1986-01-01    2016-09-12 2020_dataset           6
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                            1                             1         Blood
    3                            3                             3         Blood
    4                            7                             7         Blood
    5                           11                            14         Blood
    6                           32                            30         Blood
      visit
    1     1
    2     2
    3     3
    4     4
    5     5
    6     6

> Q10. Now using the same procedure join meta with titer data so we can
> further analyze this data in terms of time of visit aP/wP, male/female
> etc.

Let’s do one more join to join `ab_titer` with all our `meta` data.

``` r
abdata <- inner_join(ab_titer,meta)
```

    Joining with `by = join_by(specimen_id)`

``` r
head(abdata)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgE               FALSE   Total 1110.21154       2.493425
    2           1     IgE               FALSE   Total 2708.91616       2.493425
    3           1     IgG                TRUE      PT   68.56614       3.736992
    4           1     IgG                TRUE     PRN  332.12718       2.602350
    5           1     IgG                TRUE     FHA 1887.12263      34.050956
    6           1     IgE                TRUE     ACT    0.10000       1.000000
       unit lower_limit_of_detection subject_id infancy_vac biological_sex
    1 UG/ML                 2.096133          1          wP         Female
    2 IU/ML                29.170000          1          wP         Female
    3 IU/ML                 0.530000          1          wP         Female
    4 IU/ML                 6.205949          1          wP         Female
    5 IU/ML                 4.679535          1          wP         Female
    6 IU/ML                 2.816431          1          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    5 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    6 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                           -3                             0         Blood
    3                           -3                             0         Blood
    4                           -3                             0         Blood
    5                           -3                             0         Blood
    6                           -3                             0         Blood
      visit
    1     1
    2     1
    3     1
    4     1
    5     1
    6     1

> Q11. How many specimens (i.e. entries in abdata) do we have for each
> isotype?

``` r
table(abdata$isotype)
```


      IgE   IgG  IgG1  IgG2  IgG3  IgG4 
     6698  7265 11993 12000 12000 12000 

``` r
table(abdata$antigen)
```


        ACT   BETV1      DT   FELD1     FHA  FIM2/3   LOLP1     LOS Measles     OVA 
       1970    1970    6318    1970    6712    6318    1970    1970    1970    6318 
        PD1     PRN      PT     PTM   Total      TT 
       1970    6712    6712    1970     788    6318 

Lets focus on IgG isotype

``` r
igg <-abdata|>
      filter(isotype== "IgG")
head(igg)
```

      specimen_id isotype is_antigen_specific antigen        MFI MFI_normalised
    1           1     IgG                TRUE      PT   68.56614       3.736992
    2           1     IgG                TRUE     PRN  332.12718       2.602350
    3           1     IgG                TRUE     FHA 1887.12263      34.050956
    4          19     IgG                TRUE      PT   20.11607       1.096366
    5          19     IgG                TRUE     PRN  976.67419       7.652635
    6          19     IgG                TRUE     FHA   60.76626       1.096457
       unit lower_limit_of_detection subject_id infancy_vac biological_sex
    1 IU/ML                 0.530000          1          wP         Female
    2 IU/ML                 6.205949          1          wP         Female
    3 IU/ML                 4.679535          1          wP         Female
    4 IU/ML                 0.530000          3          wP         Female
    5 IU/ML                 6.205949          3          wP         Female
    6 IU/ML                 4.679535          3          wP         Female
                   ethnicity  race year_of_birth date_of_boost      dataset
    1 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    2 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    3 Not Hispanic or Latino White    1986-01-01    2016-09-12 2020_dataset
    4                Unknown White    1983-01-01    2016-10-10 2020_dataset
    5                Unknown White    1983-01-01    2016-10-10 2020_dataset
    6                Unknown White    1983-01-01    2016-10-10 2020_dataset
      actual_day_relative_to_boost planned_day_relative_to_boost specimen_type
    1                           -3                             0         Blood
    2                           -3                             0         Blood
    3                           -3                             0         Blood
    4                           -3                             0         Blood
    5                           -3                             0         Blood
    6                           -3                             0         Blood
      visit
    1     1
    2     1
    3     1
    4     1
    5     1
    6     1

> Q12. What are the different \$dataset values in abdata and what do you
> notice about the number of rows for the most “recent” dataset?

> Q13. Complete the following code to make a summary boxplot of Ab titer
> levels (MFI) for all antigens

Make a plot of `MFI_normalised` values for all `antigen` values.

``` r
ggplot(igg)+
  aes(MFI_normalised,antigen) +
  geom_boxplot()
```

![](Class18_files/figure-commonmark/unnamed-chunk-16-1.png)

The antigens “PT” “FIM2/3” and “FHA” appear to have the widest range of
values

> Q. Is there a difference for these responses between aP and wP
> individuals?

``` r
ggplot(igg)+
  aes(MFI_normalised,antigen, col= infancy_vac) +
  geom_boxplot()
```

![](Class18_files/figure-commonmark/unnamed-chunk-17-1.png)

``` r
ggplot(igg)+
  aes(MFI_normalised,antigen) +
  geom_boxplot() +
facet_wrap(~infancy_vac)
```

![](Class18_files/figure-commonmark/unnamed-chunk-18-1.png)

> Q. Is there a difference with time (i.e. before booster shot vs after
> booster shot)

``` r
ggplot(igg)+
  aes(MFI_normalised,antigen, col= infancy_vac) +
  geom_boxplot() +
  facet_wrap(~visit)
```

![](Class18_files/figure-commonmark/unnamed-chunk-19-1.png)

``` r
ab.PT.21 <- abdata |>
  filter(dataset == "2021_dataset",
        isotype == "IgG",  
        antigen == "PT")


  ggplot(ab.PT.21) +
    aes(x=planned_day_relative_to_boost,
        y=MFI_normalised,
        col=infancy_vac,
        group=subject_id) +
    
    geom_line(alpha=0.25, linewidth= 0.6)+
    geom_point(alpha= 0.4, size=1.5) +
    
    
    stat_summary(fun=mean,
                 geom="smooth",
                 se= FALSE,
                 linewidth= 2.8,
                 span=0.5,
                 aes(group= infancy_vac, col= infancy_vac))+
    
    geom_vline(xintercept=0, linetype="dashed") +
    geom_vline(xintercept=14, linetype="dashed") +
    
  labs(title="2021 dataset IgG PT",
       subtitle = "Dashed lines indicate day 0 (pre-boost) and 14 (apparent peak levels)")
```

    Warning in stat_summary(fun = mean, geom = "smooth", se = FALSE, linewidth =
    2.8, : Ignoring unknown parameters: `span`

![](Class18_files/figure-commonmark/unnamed-chunk-20-1.png)

( I confirmed the the professor all the questions he wanted us to
answer)
