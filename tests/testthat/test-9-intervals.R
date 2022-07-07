test_that("zoom in and zoom out works", {
    expect_equal(
        gintervals(1, 0, 200),
        gintervals(1, 0, 200) %>%
            gintervals.zoom_in(2) %>%
            gintervals.zoom_out(2)
    )
    expect_equal(gintervals(1, 0, 200) %>% gintervals.zoom_in(2), gintervals(1, 50, 150))
    expect_equal(gintervals(1, 50, 150) %>% gintervals.zoom_out(2), gintervals(1, 0, 200))
    expect_equal(gintervals(1, 0, 200) %>% gintervals.zoom_out(2), gintervals(1, 0, 300))
})

test_that("shift left and right works", {
    expect_equal(
        gintervals(1, 50, 200),
        gintervals(1, 50, 200) %>%
            gintervals.shift_left(20) %>%
            gintervals.shift_right(20)
    )
    expect_equal(gintervals(1, 0, 200) %>% gintervals.shift_left(20), gintervals(1, 0, 180))
    expect_equal(gintervals(1, 50, 200) %>% gintervals.shift_right(20), gintervals(1, 70, 220))
    expect_equal(gintervals(1, 0, 200) %>% gintervals.shift_right(20), gintervals(1, 20, 220))
})

test_that("extend left and right works", {
    expect_equal(
        gintervals(1, 30, 220),
        gintervals(1, 50, 200) %>%
            gintervals.extend_left(20) %>%
            gintervals.extend_right(20)
    )
    expect_equal(gintervals(1, 50, 200) %>% gintervals.extend_left(20), gintervals(1, 30, 200))
    expect_equal(gintervals(1, 50, 200) %>% gintervals.extend_right(20), gintervals(1, 50, 220))
    expect_equal(gintervals(1, 0, 200) %>% gintervals.extend_left(20), gintervals(1, 0, 200))
})
