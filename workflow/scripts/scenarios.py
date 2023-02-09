import demes

class Scenario:
    pulse_times = [
        200,
        180,
        160,
        140,
        120,
        100,
        80,
        60,
    ]
    
    def __init__(self, name: str, N_anc: int, pulses: list[list], filename: str):
        assert len(pulses[0]) == N_anc
        assert len(pulses) == len(self.pulse_times)
        self.name = name
        self.N_anc = N_anc
        self.pulses = pulses
        self.filename = filename

    def build(self):
        b = demes.Builder(
            description=self.name,
            time_units="generations",
            generation_time=1,
        )
        b.add_deme(
            "Pop0",
            description="Ancestral 1",
            epochs=[dict(end_time=0, start_size=5000)],
        )
        start = 1500
        for i in range(1, self.N_anc):
            b.add_deme(
                f"Pop{i}",
                description=f"Ancestral {i + 1}",
                ancestors=["Pop0"],
                start_time=start,
                epochs=[dict(end_time=0, start_size=5000)],
            )
            start -= 200
        b.add_deme(
            f"Pop{self.N_anc}",
            description="Admixed",
            ancestors=["Pop0"],
            start_time=210,
            epochs=[dict(end_time=0, start_size=5000)],
        )
        for t, p in zip(self.pulse_times, self.pulses):
            b.add_pulse(
                sources=[f"Pop{i}" for i in range(self.N_anc)],
                dest=f"Pop{self.N_anc}",
                proportions=p,
                time=t,
            )
        self.graph = b.resolve()
        demes.dump(self.graph, self.filename)

#==================
# Two populations
#==================

# Scenario 2A
S2A = Scenario(
    name="2A",
    N_anc=2,
    pulses=[
        [0, 0.1],
        [0, 0.1],
        [0, 0.1],
        [0, 0],
        [0, 0.1],
        [0, 0.1],
        [0, 0.1],
        [0, 0],
    ],
    filename=snakemake.output[0]
)
S2A.build()

# Scenario 2B
S2B = Scenario(
    name="2B",
    N_anc=2,
    pulses=[
        [0, 0.1],
        [0, 0.1],
        [0, 0.1],
        [0, 0],
        [0.1, 0],
        [0.1, 0],
        [0.1, 0],
        [0, 0],
    ],
    filename=snakemake.output[1]
)
S2B.build()

# Scenario 2C
S2C = Scenario(
    name="2C",
    N_anc=2,
    pulses=[
        [0, 0.2],
        [0.2, 0],
        [0, 0.2],
        [0, 0],
        [0.2, 0],
        [0, 0.2],
        [0.2, 0],
        [0, 0],
    ],
    filename=snakemake.output[2]
)
S2C.build()


#==================
# Three populations
#==================

# Scenario 3A
S3A = Scenario(
    name="3A",
    N_anc=3,
    pulses=[
        [0, 0.1, 0],
        [0, 0.1, 0],
        [0, 0.1, 0],
        [0, 0, 0],
        [0, 0, 0.1],
        [0, 0, 0.1],
        [0, 0, 0.1],
        [0, 0, 0],
    ],
    filename=snakemake.output[3]
)
S3A.build()

# Scenario 3B
S3B = Scenario(
    name="3B",
    N_anc=3,
    pulses=[
        [0, 0.2, 0],
        [0, 0, 0.2],
        [0, 0.2, 0],
        [0, 0, 0],
        [0, 0, 0.2],
        [0, 0.2, 0],
        [0, 0, 0.2],
        [0, 0, 0],
    ],
    filename=snakemake.output[4]
)
S3B.build()
