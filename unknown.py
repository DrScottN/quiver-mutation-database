# Code to handle when our properties are unknown rather than just true or false
# Goal is to have something a bit more robust than None to use


class unknown:
    def __init__(self):
        pass

    def __str__(self):
        return "Unknown"

    def __eq__(self, other):
        return isinstance(other, unknown)

    def __and__(self, other):
        if isinstance(other, unknown):
            return self

        if other:
            return self

        return other

    def __rand__(self, other):
        return self.__and__(other)

    def __or__(self, other):
        if isinstance(other, unknown):
            return self

        if not other:
            return self

        return other

    def __ror__(self,other):
        return self.__or__(other)

    def __bool__(self):
        raise Exception("You just tried to convert Unknown to a Boolean. Can't do that!")

        return self

    def __invert__(self):
        return self

def test():
    Values = [Unknown, True, False]

    for a in Values[:1]:
        for b in Values:
            print("-----------------")
            print(a, b)
            print(a & b)
            print(a | b)
            print(a == b)
            print(~b)

# Can now import unknown.Unknown as Unknown and use it as desired
Unknown = unknown()

if __name__ == "__main__":
    test()
