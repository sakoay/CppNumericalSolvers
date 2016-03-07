package(default_visibility = ["//visibility:public"])


cc_library(
    name = "cppnumericalsolvers",
    hdrs = glob(["include/cppoptlib/**/*.h"]),
    includes = [ includes ],
    visibility = ["//visibility:public"],
)


cc_binary(
    name = "linear regression",
    srcs = ["hsrc/examples/linearregression.cpp"],
    deps = [":cppnumericalsolvers"],
)
