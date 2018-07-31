def assertlength(inputset, expectedlength):
    message = "Assertion Length Error, expected: " + str(expectedlength) + " actual:" + str(len(inputset))
    #print(message)
    assert (len(inputset) == expectedlength)
