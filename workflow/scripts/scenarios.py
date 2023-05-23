import demes
import os
import demesdraw

class Scenario:
    pulse_times = [
        150,
        130,
        110,
        90,
        70,
        50,
        30,
        10,
    ]
    
    def __init__(
        self,
        name: str,
        N_anc: int,
        pop_sizes: list[int],
        pulses: list[list],
        file_prefix: str
    ):
        assert len(pulses[0]) == N_anc
        assert len(pop_sizes) == (N_anc + 1)
        assert len(pulses) == len(self.pulse_times)
        self.name = name
        self.N_anc = N_anc
        self.pop_sizes = pop_sizes
        self.pulses = pulses
        self.file_prefix = file_prefix
        self.plot = f"{self.file_prefix}.svg"

    def build(self):
        b = demes.Builder(
            description=self.name,
            time_units="generations",
            generation_time=1,
        )
        b.add_deme(
            "Pop0",
            description="Ancestral 1",
            epochs=[dict(end_time=0, start_size=pop_sizes[0])],
        )
        start = 1500
        for i in range(1, self.N_anc):
            b.add_deme(
                f"Pop{i}",
                description=f"Ancestral {i + 1}",
                ancestors=["Pop0"],
                start_time=start,
                epochs=[dict(end_time=0, start_size=pop_sizes[i])],
            )
            start -= 200
        b.add_deme(
            f"Pop{self.N_anc}",
            description="Admixed",
            ancestors=["Pop0"],
            start_time=200,
            epochs=[dict(end_time=0, start_size=pop_sizes[N_anc])],
        )
        for t, p in zip(self.pulse_times, self.pulses):
            b.add_pulse(
                sources=[f"Pop{i}" for i in range(self.N_anc)],
                dest=f"Pop{self.N_anc}",
                proportions=p,
                time=t,
            )
        self.graph = b.resolve()
        demes.dump(self.graph, self.file_prefix + '.yaml')
        demes.dump(self.graph, self.file_prefix + '.json', format='json', simplified=False)
        ax = demesdraw.tubes(self.graph, log_time=True)
        ax.figure.savefig(self.plot)


# ensure pulses 0s are floats!

# Scenario 2A
S2A = Scenario(
    name="2A",
    N_anc=2,
    pop_sizes=[10_000, 10_000, 10_000],
    pulses=[
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
        [.0, .0],
        [.0, .0],
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
    ],
    file_prefix=os.path.splitext(snakemake.output[0])[0],
)
S2A.build()

# Scenario 2B
S2A = Scenario(
    name="2B",
    N_anc=2,
    pop_sizes=[20_000, 1000, 10_000],
    pulses=[
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
        [.0, .0],
        [.0, .0],
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
    ],
    file_prefix=os.path.splitext(snakemake.output[0])[0],
)
S2A.build()

# Scenario 3A
S3A = Scenario(
    name="3A",
    N_anc=3,
    pop_sizes=[10_000, 10_000, 10_000, 10_000],
    pulses=[
        [.0, .2, .0],
        [.0, .2, .0],
        [.0, .2, .0],
        [.2, .0, .0],
        [.0, .0, .0],
        [.0, .0, .2],
        [.0, .0, .2],
        [.0, .0, .2],
    ],
    file_prefix=os.path.splitext(snakemake.output[1])[0],
)
S3A.build()

# Scenario 3B
S3B = Scenario(
    name="3B",
    N_anc=3,
    pop_sizes=[10_000, 10_000, 10_000, 10_000],
    pulses=[
        [.0, .2, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.2, .0, .0],
        [.0, .0, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.0, .0, .2],
    ],
    file_prefix=os.path.splitext(snakemake.output[1])[0],
)
S3B.build()

# Scenario 3C
S3C = Scenario(
    name="3C",
    N_anc=3,
    pop_sizes=[10_000, 1000, 5000, 10_000],
    pulses=[
        [.0, .2, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.2, .0, .0],
        [.0, .0, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.0, .0, .2],
    ],
    file_prefix=os.path.splitext(snakemake.output[1])[0],
)
S3C.build()

# # Scenario 2C
# S2C = Scenario(
#     name="2C",
#     N_anc=2,
#     pulses=[
#         [0, 0.2],
#         [0.2, 0],
#         [0, 0.2],
#         [0, 0],
#         [0.2, 0],
#         [0, 0.2],
#         [0.2, 0],
#         [0, 0],
#     ],
#     file_prefix =snakemake.output[2]
# )
# S2C.build()

# # Scenario 3B
# S3B = Scenario(
#     name="3B",
#     N_anc=3,
#     pulses=[
#         [0, 0.2, 0],
#         [0, 0, 0.2],
#         [0, 0.2, 0],
#         [0, 0, 0],
#         [0, 0, 0.2],
#         [0, 0.2, 0],
#         [0, 0, 0.2],
#         [0, 0, 0],
#     ],
#     file_prefix =snakemake.output[4]
# )
# S3B.build()
