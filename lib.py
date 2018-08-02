def assertlength(inputset, expectedlength,particle=None):
    message = "Assertion Length Error, expected: " + str(expectedlength) + " actual:" + str(len(inputset))
    try:
        assert (len(inputset) == expectedlength)
    except:
        print(message)
        if particle is not None:
            particle.summarize()
        assert (len(inputset) == expectedlength)
