library(RUnit)

testsuite.interpolate <- defineTestSuite("interpolate",
                                         dirs = file.path("test/interpolate"),
                                         testFileRegexp = "^test_.+\\.R",
                                         testFuncRegexp = "^test.+",
                                         rngKind = "Marsaglia-Multicarry",
                                         rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite.interpolate)
