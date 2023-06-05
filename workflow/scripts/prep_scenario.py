import demes
import os
import demesdraw

sc = snakemake.wildcards['sc']

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
        path: str
    ):
        assert len(pulses[0]) == N_anc
        assert len(pop_sizes) == (N_anc + 1)
        assert len(pulses) == len(self.pulse_times)
        self.name = name
        self.N_anc = N_anc
        self.pop_sizes = pop_sizes
        self.pulses = pulses
        self.file_prefix = path + '/scenario_' + name
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
            epochs=[dict(end_time=0, start_size=self.pop_sizes[0])],
        )
        start = 1500
        for i in range(1, self.N_anc):
            b.add_deme(
                f"Pop{i}",
                description=f"Ancestral {i + 1}",
                ancestors=["Pop0"],
                start_time=start,
                epochs=[dict(end_time=0, start_size=self.pop_sizes[i])],
            )
            start -= 200
        b.add_deme(
            f"Pop{self.N_anc}",
            description="Admixed",
            ancestors=["Pop0"],
            start_time=200,
            epochs=[dict(end_time=0, start_size=self.pop_sizes[self.N_anc])],
        )
        for t, p in zip(self.pulse_times, self.pulses):
            if sum(p) != 0: # only add pulse if there is one
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

sc_dict = dict()
# Scenario 2NGF (No Gene Flow)
sc_dict['2NGF'] = Scenario( 
    name="2NGF",
    N_anc=2,
    pop_sizes=[5_000, 5_000, 5_000],
    pulses=[
        [.0, .0],
        [.0, .0],
        [.0, .0],
        [.0, .0],
        [.0, .0],
        [.0, .0],
        [.0, .0],
        [.0, .0],
    ],
    path=snakemake.params['outdir'],
)

# Scenario 2A
sc_dict['2A'] = Scenario(
    name="2A",
    N_anc=2,
    pop_sizes=[5_000, 5_000, 5_000],
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
    path=snakemake.params['outdir'],
)

# Scenario 2B
sc_dict['2B'] = Scenario(
    name="2B",
    N_anc=2,
    pop_sizes=[5_000, 5_000, 5_000],
    pulses=[
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
        [.2, .0],
        [.2, .0],
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
    ],
    path=snakemake.params['outdir'],
)

# Scenario 2C
sc_dict['2C'] = Scenario(
    name="2C",
    N_anc=2,
    pop_sizes=[10_000, 1_000, 5_000],
    pulses=[
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
        [.2, .0],
        [.2, .0],
        [.0, 0.2],
        [.0, 0.2],
        [.0, 0.2],
    ],
    path=snakemake.params['outdir'],
)

# Scenario 3A
sc_dict['3A'] = Scenario(
    name="3A",
    N_anc=3,
    pop_sizes=[5_000, 5_000, 5_000, 5_000],
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
    path=snakemake.params['outdir'],
)

# Scenario 3B
sc_dict['3B'] = Scenario(
    name="3B",
    N_anc=3,
    pop_sizes=[5_000, 5_000, 5_000, 5_000],
    pulses=[
        [.0, .2, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.2, .0, .0],
        [.2, .0, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.0, .0, .2],
    ],
    path=snakemake.params['outdir'],
)

# Scenario 3C
sc_dict['3C'] = Scenario(
    name="3C",
    N_anc=3,
    pop_sizes=[5_000, 1_000, 10_000, 5_000],
    pulses=[
        [.0, .2, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.2, .0, .0],
        [.2, .0, .0],
        [.0, .0, .2],
        [.0, .2, .0],
        [.0, .0, .2],
    ],
    path=snakemake.params['outdir'],
)


sc_dict[sc].build()