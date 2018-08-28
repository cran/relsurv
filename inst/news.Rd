\name{NEWS}
\title{NEWS file for the relsurv package}

\section{Changes in version 2.2-1}{
  \itemize{
    \item 10 Aug 2018  Corrected a bug in rformulate. R in (rtable)date format is put into rform$data, the original format of the variables is not preserved
    }
}
\section{Changes in version 2.2}{
  \itemize{
    \item 15 Apr 2018  Multiple changes to rformulate function (by Terry Therneau) to be in line with the new survival package requirements - several date formats are now allowed (date, Date, POSIXt)
    \item 7 Aug 2018   Add the rmap argument to functions rs.surv, rsmul, rsadd, rstrans, nessie, rs.period, rsdiff,cmp.rel, as is the practice in the survival package, and update the manual pages and examples. The ratetable() argument in the formula is still allowed but flagged as deprecated.
    \item Allow all the transrate functions to work without the dimid attribute
    \item New Slovene population tables included (up to 2016)
    
}

}

