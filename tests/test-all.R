library(testthat)
library(phest)

test_that("Weibull estimates are sensible",{

    # These are the Roberts & Solow data (2003)
    expect_equal(
        weib.limit(c(1662, 1638, 1631, 1628, 1628, 1611, 1607, 1602, 1601, 1598), upper=TRUE),
        setNames(c(1690, 1799, 1669), c("estimate","conf-interval","conf-interval")),
        tolerance=.1
    )
    expect_equal(
        weib.limit(c(1662, 1638, 1631, 1628, 1628, 1611, 1607, 1602, 1601, 1598)),
        setNames(c(1595, 1575, 1597), c("estimate","conf-interval","conf-interval")),
        tolerance=.1
    )
    
})
