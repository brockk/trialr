
devtools::install_github("hadley/pkgdown")

pkgdown::build_site()
# Breaks:
# Error in .Call("CPP_stanc280", model_code, model_cppname, allow_undefined,  (from stanmodels.R#31) :
# "CPP_stanc280" not resolved from current namespace (rstan)
? pkgdown::build_site

# build_site wraps five functions:
# init_site()
# build_articles()
# build_home()
# build_reference()
# build_news()

# The first seems to have at-least-semi-succeeded
# Let's try the others
pkgdown::build_articles()
# This seemed to succeed.

pkgdown::build_home()
# This seemed to succeed.

pkgdown::build_reference()
# This breaks

pkgdown::build_news()
# This seemed to succeed.

? pkgdown::build_reference
