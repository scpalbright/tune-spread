
class VerbosityHandler:

    def __init__(self, level: int = 1):
        self.level = level

    def __call__(self, msg: str, level: int = 1):
        if level <= self.level:
            print(msg)
