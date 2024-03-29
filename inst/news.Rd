\name{NEWS}
\title{NEWS file for the relsurv package}

\section{Changes in version 2.2-9}{
  \itemize{
    \item 21 December 2022  Make a few improvements in the years function, primarily for measure='yl2017'
    \item 21 December 2022  Implement left-truncation in rsadd (time-dependent covariates already available)
}
}

\section{Changes in version 2.2-8}{
  \itemize{
    \item 12 August 2022  Update Sint objects for new version
}
}

\section{Changes in version 2.2-7}{
  \itemize{
    \item 09 March 2022  The slopop object has been updated containing data for 2020
    \item 09 March 2022  The C scripts cmpfast.c, netfastpinter2.c, netfastp.c have been
    updated to correctly assign individuals at risk for left truncated data
    \item 09 March 2022  A small bug was fixed in the years function
    \item 09 March 2022  Change if(class) to if(inherits) function where Debian servers gave a note
}
}
\section{Changes in version 2.2-6}{
  \itemize{
    \item 1 December 2021  Function years has been added for calculating the years life difference
}
}
\section{Changes in version 2.2-5}{
  \itemize{
    \item 5 July 2021  exp.prep is copied and exported as expprep2
    \item 5 July 2021  All R code has been roxygenized
}
}
\section{Changes in version 2.2-4}{
  \itemize{
    \item 8 June 2021  The latest version of slopop is added (1930-2019)
    \item 8 June 2021  The C code routine netfastpinter has been changed to correctly assign individuals at risk for left truncated data
    \item 8 June 2021  In the summary of ratetable objects as.Date is used (instead of as.date)
    \item 8 June 2021  A C code routine netweiDM has been added (based on netwei)
    }
}


\section{Changes in version 2.2-3}{
  \itemize{
    \item 28 Nov 2018  The CITATION changed to include the paper descrbing the package published in JSS
    }
}


\section{Changes in version 2.2-2}{
  \itemize{
    \item 10 Oct 2018  Corrected a bug in rformulate. Strata did not work correctly.
    \item 16 Oct 2018  Removed package splines from Depends to Imports. Set the depends for package survival to >= 2.42
    }
}


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

