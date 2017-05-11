library(RUnit)

testsuite.logspec <- defineTestSuite("logspec",
                                     dirs = file.path("test/logspec"),
                                     testFileRegexp = "^test_.+\\.R",
                                     testFuncRegexp = "^test.+",
                                     rngKind = "Marsaglia-Multicarry",
                                     rngNormalKind = "Kinderman-Ramage")

testResult <- runTestSuite(testsuite.logspec)

printTextProtocol(testResult, showDetails = TRUE)
