"invtime" <-
function (y = 0.1, age = 23011, sex = "male", year = 9497, scale = 1, 
    ratetable = slop, lower, upper) 
{
    if (!is.numeric(age)) 
        stop("\"age\" must be numeric", call. = FALSE)
    if (!is.numeric(y)) 
        stop("\"y\" must be numeric", call. = FALSE)
    if (!is.numeric(scale)) 
        stop("\"scale\" must be numeric", call. = FALSE)
    temp <- data.frame(age = age, sex = I(sex), year = year)
    if (missing(lower)) {
        if (!missing(upper)) 
            stop("Argument \"lower\" is missing, with no default", 
                call. = FALSE)
        nyears <- round((110 - age/365.24))
        tab <- data.frame(age = rep(age, nyears), sex = I(rep(sex, 
            nyears)), year = rep(year, nyears))
        vred <- 1 - survexp(c(0, 1:(nyears - 1)) * 365.24 ~ ratetable(age = age, 
            sex = sex, year = year), ratetable = ratetable, data = tab, 
            cohort = FALSE)
        place <- sum(vred <= y)
        if (place == 0) 
            lower <- 0
        else lower <- floor((place - 1) * 365.24 - place)
        upper <- ceiling(place * 365.24 + place)
    }
    else {
        if (missing(upper)) 
            stop("Argument \"upper\" is missing, with no default", 
                call. = FALSE)
        if (!is.integer(lower)) 
            lower <- floor(lower)
        if (!is.integer(upper)) 
            upper <- ceiling(upper)
        if (upper <= lower) 
            stop("'upper' must be higher than 'lower'", call. = FALSE)
    }
    lower <- max(0, lower)
    tab <- data.frame(age = rep(age, upper - lower + 1), sex = I(rep(sex, 
        upper - lower + 1)), year = rep(year, upper - lower + 
        1))
    vred <- 1 - survexp((lower:upper) ~ ratetable(age = age, 
        sex = sex, year = year), ratetable = ratetable, data = tab, 
        cohort = FALSE)
    place <- sum(vred <= y)
    if (place == 0) 
        warning(paste("The event happened on or before day", 
            lower), call. = FALSE)
    if (place == length(vred)) 
        warning(paste("The event happened on or after day", upper), 
            call. = FALSE)
    t <- (place + lower-1)/scale
    age <- round(age/365.24, 0.01)
    return(list(age, sex, year, Y = y, T = t))
}
