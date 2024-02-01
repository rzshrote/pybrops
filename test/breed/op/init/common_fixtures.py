from pybrops.breed.op.init.InitializationOperator import InitializationOperator


class DummyInitializationOperator(InitializationOperator):
    def initialize(self, miscout: dict | None, **kwargs: dict) -> tuple:
        return super().initialize(miscout, **kwargs)

