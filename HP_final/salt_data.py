from dataclasses import dataclass

@dataclass
class SaltData:
    salt: str
    molality: float
    molarity: float
    num_particles: float
    osmotic_coefficient: float
    density: str

salt_infos = [
    {
        "salt": "CsBr",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.917,
        "density": "1.01356"
    },
    {
        "salt": "CsBr",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.896,
        "density": "1.02986"
    },
    {
        "salt": "CsBr",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.882,
        "density": "1.04598"
    },
    {
        "salt": "CsBr",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.873,
        "density": "1.06192"
    },
    {
        "salt": "CsBr",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.865,
        "density": "1.07769"
    },
    {
        "salt": "CsBr",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.861,
        "density": "1.09329"
    },
    {
        "salt": "CsBr",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.857,
        "density": "1.10873"
    },
    {
        "salt": "CsBr",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.854,
        "density": "1.12401"
    },
    {
        "salt": "CsBr",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.852,
        "density": "1.13912"
    },
    {
        "salt": "CsBr",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.85,
        "density": "1.15409"
    },
    {
        "salt": "CsBr",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.849,
        "density": "1.18358"
    },
    {
        "salt": "CsBr",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.848,
        "density": "1.21243"
    },
    {
        "salt": "CsBr",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.848,
        "density": "1.24073"
    },
    {
        "salt": "CsBr",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.85,
        "density": "1.26847"
    },
    {
        "salt": "CsBr",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.852,
        "density": "1.29566"
    },
    {
        "salt": "CsBr",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.859,
        "density": "1.36136"
    },
    {
        "salt": "CsBr",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.866,
        "density": "1.42396"
    },
    {
        "salt": "CsBr",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.874,
        "density": "1.48367"
    },
    {
        "salt": "CsBr",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.884,
        "density": "1.54067"
    },
    {
        "salt": "CsBr",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.892,
        "density": "1.59514"
    },
    {
        "salt": "CsBr",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.901,
        "density": "1.64724"
    },
    {
        "salt": "CsCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.917,
        "density": "1.00983"
    },
    {
        "salt": "CsCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.897,
        "density": "1.02246"
    },
    {
        "salt": "CsCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.885,
        "density": "1.03496"
    },
    {
        "salt": "CsCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.875,
        "density": "1.04732"
    },
    {
        "salt": "CsCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.869,
        "density": "1.05956"
    },
    {
        "salt": "CsCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.864,
        "density": "1.07168"
    },
    {
        "salt": "CsCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.861,
        "density": "1.08368"
    },
    {
        "salt": "CsCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.859,
        "density": "1.09557"
    },
    {
        "salt": "CsCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.858,
        "density": "1.10734"
    },
    {
        "salt": "CsCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.857,
        "density": "1.119"
    },
    {
        "salt": "CsCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.856,
        "density": "1.14201"
    },
    {
        "salt": "CsCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.856,
        "density": "1.16456"
    },
    {
        "salt": "CsCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.857,
        "density": "1.18671"
    },
    {
        "salt": "CsCl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.859,
        "density": "1.20845"
    },
    {
        "salt": "CsCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.864,
        "density": "1.22979"
    },
    {
        "salt": "CsCl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.871,
        "density": "1.28149"
    },
    {
        "salt": "CsCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.88,
        "density": "1.33093"
    },
    {
        "salt": "CsCl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.891,
        "density": "1.37823"
    },
    {
        "salt": "CsCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.901,
        "density": "1.42352"
    },
    {
        "salt": "CsCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.913,
        "density": "1.46694"
    },
    {
        "salt": "CsCl",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.923,
        "density": "1.50858"
    },
    {
        "salt": "CsCl",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 0.934,
        "density": "1.54855"
    },
    {
        "salt": "CsCl",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 0.945,
        "density": "1.58697"
    },
    {
        "salt": "CsI",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.916,
        "density": "1.01706"
    },
    {
        "salt": "CsI",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.895,
        "density": "1.0368"
    },
    {
        "salt": "CsI",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.88,
        "density": "1.05627"
    },
    {
        "salt": "CsI",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.87,
        "density": "1.0755"
    },
    {
        "salt": "CsI",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.863,
        "density": "1.09449"
    },
    {
        "salt": "CsI",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.858,
        "density": "1.11325"
    },
    {
        "salt": "CsI",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.855,
        "density": "1.13177"
    },
    {
        "salt": "CsI",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.852,
        "density": "1.15008"
    },
    {
        "salt": "CsI",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.849,
        "density": "1.16816"
    },
    {
        "salt": "CsI",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.846,
        "density": "1.18602"
    },
    {
        "salt": "CsI",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.842,
        "density": "1.22114"
    },
    {
        "salt": "CsI",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.839,
        "density": "1.25539"
    },
    {
        "salt": "CsI",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.836,
        "density": "1.28886"
    },
    {
        "salt": "CsI",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.834,
        "density": "1.32155"
    },
    {
        "salt": "CsI",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.832,
        "density": "1.35349"
    },
    {
        "salt": "CsI",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.827,
        "density": "1.43016"
    },
    {
        "salt": "CsI",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.822,
        "density": "1.50251"
    },
    {
        "salt": "CsNO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.902,
        "density": "1.01131"
    },
    {
        "salt": "CsNO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.869,
        "density": "1.02536"
    },
    {
        "salt": "CsNO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.842,
        "density": "1.03922"
    },
    {
        "salt": "CsNO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.82,
        "density": "1.0529"
    },
    {
        "salt": "CsNO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.802,
        "density": "1.06642"
    },
    {
        "salt": "CsNO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.787,
        "density": "1.07978"
    },
    {
        "salt": "CsNO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.774,
        "density": "1.09298"
    },
    {
        "salt": "CsNO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.761,
        "density": "1.10604"
    },
    {
        "salt": "CsNO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.748,
        "density": "1.11895"
    },
    {
        "salt": "CsNO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.736,
        "density": "1.13173"
    },
    {
        "salt": "CsNO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.715,
        "density": "1.15691"
    },
    {
        "salt": "CsNO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.695,
        "density": "1.18157"
    },
    {
        "salt": "KBr",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.928,
        "density": "1.00548"
    },
    {
        "salt": "KBr",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.916,
        "density": "1.0138"
    },
    {
        "salt": "KBr",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.91,
        "density": "1.02203"
    },
    {
        "salt": "KBr",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.906,
        "density": "1.03018"
    },
    {
        "salt": "KBr",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.904,
        "density": "1.03826"
    },
    {
        "salt": "KBr",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.904,
        "density": "1.04625"
    },
    {
        "salt": "KBr",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.904,
        "density": "1.05418"
    },
    {
        "salt": "KBr",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.905,
        "density": "1.06203"
    },
    {
        "salt": "KBr",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.906,
        "density": "1.06981"
    },
    {
        "salt": "KBr",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.907,
        "density": "1.07752"
    },
    {
        "salt": "KBr",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.91,
        "density": "1.09275"
    },
    {
        "salt": "KBr",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.914,
        "density": "1.10767"
    },
    {
        "salt": "KBr",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.917,
        "density": "1.12236"
    },
    {
        "salt": "KBr",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.922,
        "density": "1.1368"
    },
    {
        "salt": "KBr",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.927,
        "density": "1.15098"
    },
    {
        "salt": "KBr",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.941,
        "density": "1.18541"
    },
    {
        "salt": "KBr",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.955,
        "density": "1.21841"
    },
    {
        "salt": "KBr",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.969,
        "density": "1.25008"
    },
    {
        "salt": "KBr",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.984,
        "density": "1.28047"
    },
    {
        "salt": "KBr",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.0,
        "density": "1.30968"
    },
    {
        "salt": "KBr",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 1.015,
        "density": "1.33776"
    },
    {
        "salt": "KBr",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 1.028,
        "density": "1.36478"
    },
    {
        "salt": "KCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.927,
        "density": "1.00174"
    },
    {
        "salt": "KCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.913,
        "density": "1.00636"
    },
    {
        "salt": "KCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.906,
        "density": "1.01091"
    },
    {
        "salt": "KCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.902,
        "density": "1.01541"
    },
    {
        "salt": "KCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.899,
        "density": "1.01987"
    },
    {
        "salt": "KCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.898,
        "density": "1.02427"
    },
    {
        "salt": "KCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.897,
        "density": "1.02863"
    },
    {
        "salt": "KCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.897,
        "density": "1.03294"
    },
    {
        "salt": "KCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.897,
        "density": "1.03722"
    },
    {
        "salt": "KCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.897,
        "density": "1.04145"
    },
    {
        "salt": "KCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.899,
        "density": "1.04981"
    },
    {
        "salt": "KCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.901,
        "density": "1.05797"
    },
    {
        "salt": "KCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.904,
        "density": "1.06601"
    },
    {
        "salt": "KCl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.908,
        "density": "1.0739"
    },
    {
        "salt": "KCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.912,
        "density": "1.08166"
    },
    {
        "salt": "KCl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.924,
        "density": "1.10045"
    },
    {
        "salt": "KCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.937,
        "density": "1.11846"
    },
    {
        "salt": "KCl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.95,
        "density": "1.13572"
    },
    {
        "salt": "KCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.965,
        "density": "1.15229"
    },
    {
        "salt": "KCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.98,
        "density": "1.16822"
    },
    {
        "salt": "KI",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.932,
        "density": "1.009"
    },
    {
        "salt": "KI",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.922,
        "density": "1.02081"
    },
    {
        "salt": "KI",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.918,
        "density": "1.03248"
    },
    {
        "salt": "KI",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.917,
        "density": "1.04403"
    },
    {
        "salt": "KI",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.917,
        "density": "1.05545"
    },
    {
        "salt": "KI",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.918,
        "density": "1.06676"
    },
    {
        "salt": "KI",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.919,
        "density": "1.07795"
    },
    {
        "salt": "KI",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.922,
        "density": "1.08902"
    },
    {
        "salt": "KI",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.924,
        "density": "1.09998"
    },
    {
        "salt": "KI",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.926,
        "density": "1.11082"
    },
    {
        "salt": "KI",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.931,
        "density": "1.1322"
    },
    {
        "salt": "KI",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.937,
        "density": "1.15312"
    },
    {
        "salt": "KI",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.943,
        "density": "1.17363"
    },
    {
        "salt": "KI",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.95,
        "density": "1.19374"
    },
    {
        "salt": "KI",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.957,
        "density": "1.21346"
    },
    {
        "salt": "KI",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.974,
        "density": "1.26111"
    },
    {
        "salt": "KI",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.99,
        "density": "1.30655"
    },
    {
        "salt": "KI",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.006,
        "density": "1.34992"
    },
    {
        "salt": "KI",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.021,
        "density": "1.39139"
    },
    {
        "salt": "KI",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.032,
        "density": "1.43108"
    },
    {
        "salt": "KNO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.906,
        "density": "1.00323"
    },
    {
        "salt": "KNO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.873,
        "density": "1.00932"
    },
    {
        "salt": "KNO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.851,
        "density": "1.01532"
    },
    {
        "salt": "KNO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.833,
        "density": "1.02124"
    },
    {
        "salt": "KNO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.817,
        "density": "1.02709"
    },
    {
        "salt": "KNO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.802,
        "density": "1.03286"
    },
    {
        "salt": "KNO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.79,
        "density": "1.03857"
    },
    {
        "salt": "KNO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.778,
        "density": "1.04421"
    },
    {
        "salt": "KNO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.767,
        "density": "1.04979"
    },
    {
        "salt": "KNO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.756,
        "density": "1.0553"
    },
    {
        "salt": "KNO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.736,
        "density": "1.06616"
    },
    {
        "salt": "KNO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.718,
        "density": "1.07675"
    },
    {
        "salt": "KNO3",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.7,
        "density": "1.08714"
    },
    {
        "salt": "KNO3",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.684,
        "density": "1.09731"
    },
    {
        "salt": "KNO3",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.669,
        "density": "1.10728"
    },
    {
        "salt": "KNO3",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.631,
        "density": "1.13135"
    },
    {
        "salt": "KNO3",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.602,
        "density": "1.15432"
    },
    {
        "salt": "KNO3",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.577,
        "density": "1.1763"
    },
    {
        "salt": "LiBr",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.943,
        "density": "1.00327"
    },
    {
        "salt": "LiBr",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.944,
        "density": "1.00942"
    },
    {
        "salt": "LiBr",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.952,
        "density": "1.01552"
    },
    {
        "salt": "LiBr",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.96,
        "density": "1.02156"
    },
    {
        "salt": "LiBr",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.97,
        "density": "1.02756"
    },
    {
        "salt": "LiBr",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.981,
        "density": "1.03352"
    },
    {
        "salt": "LiBr",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.993,
        "density": "1.03944"
    },
    {
        "salt": "LiBr",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 1.007,
        "density": "1.04532"
    },
    {
        "salt": "LiBr",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 1.021,
        "density": "1.05116"
    },
    {
        "salt": "LiBr",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 1.035,
        "density": "1.05697"
    },
    {
        "salt": "LiBr",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.067,
        "density": "1.06848"
    },
    {
        "salt": "LiBr",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.098,
        "density": "1.07982"
    },
    {
        "salt": "LiBr",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.13,
        "density": "1.09104"
    },
    {
        "salt": "LiBr",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.163,
        "density": "1.10213"
    },
    {
        "salt": "LiBr",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.196,
        "density": "1.11309"
    },
    {
        "salt": "LiBr",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.278,
        "density": "1.13996"
    },
    {
        "salt": "LiBr",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.364,
        "density": "1.16608"
    },
    {
        "salt": "LiBr",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.467,
        "density": "1.1915"
    },
    {
        "salt": "LiBr",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.578,
        "density": "1.21626"
    },
    {
        "salt": "LiBr",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.687,
        "density": "1.24037"
    },
    {
        "salt": "LiBr",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 1.793,
        "density": "1.26389"
    },
    {
        "salt": "LiBr",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 1.891,
        "density": "1.28682"
    },
    {
        "salt": "LiBr",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 1.989,
        "density": "1.3092"
    },
    {
        "salt": "LiCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.939,
        "density": "0.99953"
    },
    {
        "salt": "LiCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.939,
        "density": "1.00197"
    },
    {
        "salt": "LiCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.945,
        "density": "1.00437"
    },
    {
        "salt": "LiCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.954,
        "density": "1.00675"
    },
    {
        "salt": "LiCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.963,
        "density": "1.0091"
    },
    {
        "salt": "LiCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.973,
        "density": "1.01143"
    },
    {
        "salt": "LiCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.984,
        "density": "1.01373"
    },
    {
        "salt": "LiCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.995,
        "density": "1.01602"
    },
    {
        "salt": "LiCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 1.006,
        "density": "1.01828"
    },
    {
        "salt": "LiCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 1.018,
        "density": "1.02052"
    },
    {
        "salt": "LiCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.041,
        "density": "1.02497"
    },
    {
        "salt": "LiCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.066,
        "density": "1.02931"
    },
    {
        "salt": "LiCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.091,
        "density": "1.0336"
    },
    {
        "salt": "LiCl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.116,
        "density": "1.03782"
    },
    {
        "salt": "LiCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.142,
        "density": "1.04198"
    },
    {
        "salt": "LiCl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.212,
        "density": "1.05213"
    },
    {
        "salt": "LiCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.286,
        "density": "1.06193"
    },
    {
        "salt": "LiCl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.366,
        "density": "1.07141"
    },
    {
        "salt": "LiCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.449,
        "density": "1.08061"
    },
    {
        "salt": "LiCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.533,
        "density": "1.08956"
    },
    {
        "salt": "LiCl",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 1.619,
        "density": "1.09828"
    },
    {
        "salt": "LiCl",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 1.705,
        "density": "1.10679"
    },
    {
        "salt": "LiCl",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 1.791,
        "density": "1.11511"
    },
    {
        "salt": "LiClO4",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.951,
        "density": "1.00329"
    },
    {
        "salt": "LiClO4",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.959,
        "density": "1.00944"
    },
    {
        "salt": "LiClO4",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.971,
        "density": "1.01552"
    },
    {
        "salt": "LiClO4",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.985,
        "density": "1.02153"
    },
    {
        "salt": "LiClO4",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.999,
        "density": "1.02748"
    },
    {
        "salt": "LiClO4",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 1.013,
        "density": "1.03337"
    },
    {
        "salt": "LiClO4",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 1.027,
        "density": "1.0392"
    },
    {
        "salt": "LiClO4",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 1.043,
        "density": "1.04497"
    },
    {
        "salt": "LiClO4",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 1.058,
        "density": "1.05068"
    },
    {
        "salt": "LiClO4",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 1.072,
        "density": "1.05634"
    },
    {
        "salt": "LiClO4",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.104,
        "density": "1.0675"
    },
    {
        "salt": "LiClO4",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.137,
        "density": "1.07841"
    },
    {
        "salt": "LiClO4",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.17,
        "density": "1.08913"
    },
    {
        "salt": "LiClO4",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.204,
        "density": "1.09965"
    },
    {
        "salt": "LiClO4",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.238,
        "density": "1.10998"
    },
    {
        "salt": "LiClO4",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.328,
        "density": "1.13495"
    },
    {
        "salt": "LiClO4",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.419,
        "density": "1.1588"
    },
    {
        "salt": "LiClO4",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.512,
        "density": "1.18158"
    },
    {
        "salt": "LiClO4",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.595,
        "density": "1.20337"
    },
    {
        "salt": "LiI",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.952,
        "density": "1.0068"
    },
    {
        "salt": "LiI",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.966,
        "density": "1.01645"
    },
    {
        "salt": "LiI",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.98,
        "density": "1.02602"
    },
    {
        "salt": "LiI",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.995,
        "density": "1.03552"
    },
    {
        "salt": "LiI",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 1.008,
        "density": "1.04493"
    },
    {
        "salt": "LiI",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 1.022,
        "density": "1.05427"
    },
    {
        "salt": "LiI",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 1.034,
        "density": "1.06354"
    },
    {
        "salt": "LiI",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 1.049,
        "density": "1.07274"
    },
    {
        "salt": "LiI",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 1.063,
        "density": "1.08186"
    },
    {
        "salt": "LiI",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 1.08,
        "density": "1.09092"
    },
    {
        "salt": "LiI",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.111,
        "density": "1.10885"
    },
    {
        "salt": "LiI",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.143,
        "density": "1.12649"
    },
    {
        "salt": "LiI",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.176,
        "density": "1.14388"
    },
    {
        "salt": "LiI",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.212,
        "density": "1.16103"
    },
    {
        "salt": "LiI",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.25,
        "density": "1.17794"
    },
    {
        "salt": "LiI",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.351,
        "density": "1.21922"
    },
    {
        "salt": "LiI",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.467,
        "density": "1.25914"
    },
    {
        "salt": "LiNO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.938,
        "density": "1.00104"
    },
    {
        "salt": "LiNO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.935,
        "density": "1.00499"
    },
    {
        "salt": "LiNO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.94,
        "density": "1.00889"
    },
    {
        "salt": "LiNO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.946,
        "density": "1.01277"
    },
    {
        "salt": "LiNO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.954,
        "density": "1.01661"
    },
    {
        "salt": "LiNO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.962,
        "density": "1.02041"
    },
    {
        "salt": "LiNO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.97,
        "density": "1.02419"
    },
    {
        "salt": "LiNO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.978,
        "density": "1.02794"
    },
    {
        "salt": "LiNO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.987,
        "density": "1.03167"
    },
    {
        "salt": "LiNO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.997,
        "density": "1.03536"
    },
    {
        "salt": "LiNO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.015,
        "density": "1.04268"
    },
    {
        "salt": "LiNO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.033,
        "density": "1.04986"
    },
    {
        "salt": "LiNO3",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.052,
        "density": "1.05695"
    },
    {
        "salt": "LiNO3",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.07,
        "density": "1.06394"
    },
    {
        "salt": "LiNO3",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.088,
        "density": "1.07083"
    },
    {
        "salt": "LiNO3",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.134,
        "density": "1.08764"
    },
    {
        "salt": "LiNO3",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.181,
        "density": "1.10386"
    },
    {
        "salt": "LiNO3",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.227,
        "density": "1.11953"
    },
    {
        "salt": "LiNO3",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.27,
        "density": "1.13467"
    },
    {
        "salt": "LiNO3",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.312,
        "density": "1.14931"
    },
    {
        "salt": "LiNO3",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 1.352,
        "density": "1.16347"
    },
    {
        "salt": "LiNO3",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 1.387,
        "density": "1.17717"
    },
    {
        "salt": "LiNO3",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 1.42,
        "density": "1.19043"
    },
    {
        "salt": "NaBr",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.934,
        "density": "1.00489"
    },
    {
        "salt": "NaBr",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.928,
        "density": "1.01266"
    },
    {
        "salt": "NaBr",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.928,
        "density": "1.02035"
    },
    {
        "salt": "NaBr",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.929,
        "density": "1.02799"
    },
    {
        "salt": "NaBr",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.933,
        "density": "1.03556"
    },
    {
        "salt": "NaBr",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.937,
        "density": "1.04308"
    },
    {
        "salt": "NaBr",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.942,
        "density": "1.05055"
    },
    {
        "salt": "NaBr",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.947,
        "density": "1.05796"
    },
    {
        "salt": "NaBr",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.953,
        "density": "1.06532"
    },
    {
        "salt": "NaBr",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.958,
        "density": "1.07263"
    },
    {
        "salt": "NaBr",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.969,
        "density": "1.08712"
    },
    {
        "salt": "NaBr",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.983,
        "density": "1.10138"
    },
    {
        "salt": "NaBr",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.997,
        "density": "1.11547"
    },
    {
        "salt": "NaBr",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.012,
        "density": "1.12937"
    },
    {
        "salt": "NaBr",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.028,
        "density": "1.14309"
    },
    {
        "salt": "NaBr",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.067,
        "density": "1.17662"
    },
    {
        "salt": "NaBr",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.107,
        "density": "1.2091"
    },
    {
        "salt": "NaBr",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.15,
        "density": "1.24059"
    },
    {
        "salt": "NaBr",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.199,
        "density": "1.27113"
    },
    {
        "salt": "NaCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.932,
        "density": "1.00116"
    },
    {
        "salt": "NaCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.925,
        "density": "1.0052"
    },
    {
        "salt": "NaCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.922,
        "density": "1.00921"
    },
    {
        "salt": "NaCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.92,
        "density": "1.01317"
    },
    {
        "salt": "NaCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.921,
        "density": "1.01709"
    },
    {
        "salt": "NaCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.923,
        "density": "1.02098"
    },
    {
        "salt": "NaCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.926,
        "density": "1.02483"
    },
    {
        "salt": "NaCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.929,
        "density": "1.02866"
    },
    {
        "salt": "NaCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.932,
        "density": "1.03245"
    },
    {
        "salt": "NaCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.936,
        "density": "1.03621"
    },
    {
        "salt": "NaCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.943,
        "density": "1.04366"
    },
    {
        "salt": "NaCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.951,
        "density": "1.05096"
    },
    {
        "salt": "NaCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.962,
        "density": "1.05817"
    },
    {
        "salt": "NaCl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.972,
        "density": "1.06527"
    },
    {
        "salt": "NaCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.983,
        "density": "1.07227"
    },
    {
        "salt": "NaCl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.013,
        "density": "1.08932"
    },
    {
        "salt": "NaCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.045,
        "density": "1.10579"
    },
    {
        "salt": "NaCl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.08,
        "density": "1.1217"
    },
    {
        "salt": "NaCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 1.116,
        "density": "1.13709"
    },
    {
        "salt": "NaCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 1.153,
        "density": "1.15199"
    },
    {
        "salt": "NaCl",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 1.192,
        "density": "1.16644"
    },
    {
        "salt": "NaCl",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 1.231,
        "density": "1.18048"
    },
    {
        "salt": "NaCl",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 1.271,
        "density": "1.19412"
    },
    {
        "salt": "NaI",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.938,
        "density": "1.00842"
    },
    {
        "salt": "NaI",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.936,
        "density": "1.01968"
    },
    {
        "salt": "NaI",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.939,
        "density": "1.03084"
    },
    {
        "salt": "NaI",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.945,
        "density": "1.0419"
    },
    {
        "salt": "NaI",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.952,
        "density": "1.05287"
    },
    {
        "salt": "NaI",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.959,
        "density": "1.06374"
    },
    {
        "salt": "NaI",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.967,
        "density": "1.07452"
    },
    {
        "salt": "NaI",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.975,
        "density": "1.08522"
    },
    {
        "salt": "NaI",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.983,
        "density": "1.09582"
    },
    {
        "salt": "NaI",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.991,
        "density": "1.10634"
    },
    {
        "salt": "NaI",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.007,
        "density": "1.12713"
    },
    {
        "salt": "NaI",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.025,
        "density": "1.14756"
    },
    {
        "salt": "NaI",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.043,
        "density": "1.16767"
    },
    {
        "salt": "NaI",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 1.061,
        "density": "1.18747"
    },
    {
        "salt": "NaI",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.079,
        "density": "1.20696"
    },
    {
        "salt": "NaI",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 1.129,
        "density": "1.25438"
    },
    {
        "salt": "NaI",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 1.188,
        "density": "1.30002"
    },
    {
        "salt": "NaI",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 1.243,
        "density": "1.34399"
    },
    {
        "salt": "NaNO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.921,
        "density": "1.00264"
    },
    {
        "salt": "NaNO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.902,
        "density": "1.00815"
    },
    {
        "salt": "NaNO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.89,
        "density": "1.01359"
    },
    {
        "salt": "NaNO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.881,
        "density": "1.01896"
    },
    {
        "salt": "NaNO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.873,
        "density": "1.02428"
    },
    {
        "salt": "NaNO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.867,
        "density": "1.02955"
    },
    {
        "salt": "NaNO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.862,
        "density": "1.03476"
    },
    {
        "salt": "NaNO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.858,
        "density": "1.03992"
    },
    {
        "salt": "NaNO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.854,
        "density": "1.04504"
    },
    {
        "salt": "NaNO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.851,
        "density": "1.05011"
    },
    {
        "salt": "NaNO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.845,
        "density": "1.06013"
    },
    {
        "salt": "NaNO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.839,
        "density": "1.06994"
    },
    {
        "salt": "NaNO3",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.835,
        "density": "1.0796"
    },
    {
        "salt": "NaNO3",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.83,
        "density": "1.08909"
    },
    {
        "salt": "NaNO3",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.826,
        "density": "1.09842"
    },
    {
        "salt": "NaNO3",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.817,
        "density": "1.12106"
    },
    {
        "salt": "NaNO3",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.81,
        "density": "1.14276"
    },
    {
        "salt": "NaNO3",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.804,
        "density": "1.16356"
    },
    {
        "salt": "NaNO3",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.797,
        "density": "1.1835"
    },
    {
        "salt": "NaNO3",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.792,
        "density": "1.20261"
    },
    {
        "salt": "NaNO3",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.788,
        "density": "1.22092"
    },
    {
        "salt": "NaNO3",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 0.787,
        "density": "1.23845"
    },
    {
        "salt": "NaNO3",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 0.788,
        "density": "1.25525"
    },
    {
        "salt": "RbBr",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.922,
        "density": "1.00957"
    },
    {
        "salt": "RbBr",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.904,
        "density": "1.02196"
    },
    {
        "salt": "RbBr",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.896,
        "density": "1.03421"
    },
    {
        "salt": "RbBr",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.89,
        "density": "1.04634"
    },
    {
        "salt": "RbBr",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.886,
        "density": "1.05836"
    },
    {
        "salt": "RbBr",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.884,
        "density": "1.07025"
    },
    {
        "salt": "RbBr",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.881,
        "density": "1.08203"
    },
    {
        "salt": "RbBr",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.88,
        "density": "1.0937"
    },
    {
        "salt": "RbBr",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.879,
        "density": "1.10526"
    },
    {
        "salt": "RbBr",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.878,
        "density": "1.11671"
    },
    {
        "salt": "RbBr",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.878,
        "density": "1.13931"
    },
    {
        "salt": "RbBr",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.878,
        "density": "1.16145"
    },
    {
        "salt": "RbBr",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.88,
        "density": "1.1832"
    },
    {
        "salt": "RbBr",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.882,
        "density": "1.20455"
    },
    {
        "salt": "RbBr",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.886,
        "density": "1.22551"
    },
    {
        "salt": "RbBr",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.893,
        "density": "1.2763"
    },
    {
        "salt": "RbBr",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.901,
        "density": "1.32489"
    },
    {
        "salt": "RbBr",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.911,
        "density": "1.37145"
    },
    {
        "salt": "RbBr",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.921,
        "density": "1.41611"
    },
    {
        "salt": "RbBr",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.931,
        "density": "1.45903"
    },
    {
        "salt": "RbBr",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.94,
        "density": "1.50035"
    },
    {
        "salt": "RbCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.923,
        "density": "1.00585"
    },
    {
        "salt": "RbCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.907,
        "density": "1.01453"
    },
    {
        "salt": "RbCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.898,
        "density": "1.02312"
    },
    {
        "salt": "RbCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.893,
        "density": "1.03163"
    },
    {
        "salt": "RbCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.889,
        "density": "1.04005"
    },
    {
        "salt": "RbCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.887,
        "density": "1.0484"
    },
    {
        "salt": "RbCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.886,
        "density": "1.05667"
    },
    {
        "salt": "RbCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.886,
        "density": "1.06486"
    },
    {
        "salt": "RbCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.885,
        "density": "1.07298"
    },
    {
        "salt": "RbCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.885,
        "density": "1.08104"
    },
    {
        "salt": "RbCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.886,
        "density": "1.09695"
    },
    {
        "salt": "RbCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.888,
        "density": "1.11255"
    },
    {
        "salt": "RbCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.89,
        "density": "1.12791"
    },
    {
        "salt": "RbCl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.893,
        "density": "1.14301"
    },
    {
        "salt": "RbCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.896,
        "density": "1.15785"
    },
    {
        "salt": "RbCl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.905,
        "density": "1.19391"
    },
    {
        "salt": "RbCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.916,
        "density": "1.22852"
    },
    {
        "salt": "RbCl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.928,
        "density": "1.26176"
    },
    {
        "salt": "RbCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.941,
        "density": "1.29369"
    },
    {
        "salt": "RbCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.952,
        "density": "1.32438"
    },
    {
        "salt": "RbCl",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.966,
        "density": "1.35391"
    },
    {
        "salt": "RbNO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.903,
        "density": "1.00735"
    },
    {
        "salt": "RbNO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.871,
        "density": "1.0175"
    },
    {
        "salt": "RbNO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.847,
        "density": "1.02751"
    },
    {
        "salt": "RbNO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.826,
        "density": "1.03741"
    },
    {
        "salt": "RbNO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.809,
        "density": "1.04719"
    },
    {
        "salt": "RbNO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.794,
        "density": "1.05685"
    },
    {
        "salt": "RbNO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.781,
        "density": "1.06641"
    },
    {
        "salt": "RbNO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.768,
        "density": "1.07585"
    },
    {
        "salt": "RbNO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.756,
        "density": "1.0852"
    },
    {
        "salt": "RbNO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.745,
        "density": "1.09444"
    },
    {
        "salt": "RbNO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.725,
        "density": "1.11265"
    },
    {
        "salt": "RbNO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.706,
        "density": "1.13044"
    },
    {
        "salt": "RbNO3",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.689,
        "density": "1.14788"
    },
    {
        "salt": "RbNO3",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.673,
        "density": "1.16498"
    },
    {
        "salt": "RbNO3",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.656,
        "density": "1.18174"
    },
    {
        "salt": "RbNO3",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.62,
        "density": "1.22227"
    },
    {
        "salt": "RbNO3",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.588,
        "density": "1.26101"
    },
    {
        "salt": "RbNO3",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.561,
        "density": "1.29814"
    },
    {
        "salt": "RbNO3",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.538,
        "density": "1.33386"
    },
    {
        "salt": "RbNO3",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.516,
        "density": "1.36835"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.927,
        "density": "0.99875"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.913,
        "density": "1.00041"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.906,
        "density": "1.00202"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.901,
        "density": "1.0036"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.899,
        "density": "1.00515"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.897,
        "density": "1.00668"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.896,
        "density": "1.00818"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.896,
        "density": "1.00966"
    },
    {
        "salt": "NH4Cl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.896,
        "density": "1.01111"
    },
    {
        "salt": "NH4Cl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.897,
        "density": "1.01255"
    },
    {
        "salt": "NH4Cl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.898,
        "density": "1.01538"
    },
    {
        "salt": "NH4Cl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.9,
        "density": "1.0181"
    },
    {
        "salt": "NH4Cl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.903,
        "density": "1.02077"
    },
    {
        "salt": "NH4Cl",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.906,
        "density": "1.02337"
    },
    {
        "salt": "NH4Cl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.909,
        "density": "1.02592"
    },
    {
        "salt": "NH4Cl",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.918,
        "density": "1.03201"
    },
    {
        "salt": "NH4Cl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.926,
        "density": "1.03776"
    },
    {
        "salt": "NH4Cl",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.936,
        "density": "1.04319"
    },
    {
        "salt": "NH4Cl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.945,
        "density": "1.04834"
    },
    {
        "salt": "NH4Cl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.953,
        "density": "1.05322"
    },
    {
        "salt": "NH4Cl",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.958,
        "density": "1.05786"
    },
    {
        "salt": "NH4Cl",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 0.963,
        "density": "1.06227"
    },
    {
        "salt": "NH4Cl",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 0.969,
        "density": "1.06646"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.911,
        "density": "1.00026"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.89,
        "density": "1.00341"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.876,
        "density": "1.00652"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.864,
        "density": "1.00958"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.855,
        "density": "1.0126"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.847,
        "density": "1.01559"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.84,
        "density": "1.01854"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.834,
        "density": "1.02145"
    },
    {
        "salt": "NH4NO3",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.829,
        "density": "1.02432"
    },
    {
        "salt": "NH4NO3",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.823,
        "density": "1.02716"
    },
    {
        "salt": "NH4NO3",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.813,
        "density": "1.03277"
    },
    {
        "salt": "NH4NO3",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.803,
        "density": "1.0382"
    },
    {
        "salt": "NH4NO3",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.793,
        "density": "1.04353"
    },
    {
        "salt": "NH4NO3",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.785,
        "density": "1.04873"
    },
    {
        "salt": "NH4NO3",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.776,
        "density": "1.05383"
    },
    {
        "salt": "NH4NO3",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.758,
        "density": "1.06607"
    },
    {
        "salt": "NH4NO3",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.743,
        "density": "1.07769"
    },
    {
        "salt": "NH4NO3",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.728,
        "density": "1.08871"
    },
    {
        "salt": "NH4NO3",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.715,
        "density": "1.09921"
    },
    {
        "salt": "NH4NO3",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 0.702,
        "density": "1.10921"
    },
    {
        "salt": "NH4NO3",
        "molality": 5.0,
        "molarity": 5.004,
        "num_particles": 334,
        "osmotic_coefficient": 0.69,
        "density": "1.11877"
    },
    {
        "salt": "NH4NO3",
        "molality": 5.5,
        "molarity": 5.505,
        "num_particles": 367,
        "osmotic_coefficient": 0.679,
        "density": "1.12792"
    },
    {
        "salt": "NH4NO3",
        "molality": 6.0,
        "molarity": 6.005,
        "num_particles": 400,
        "osmotic_coefficient": 0.67,
        "density": "1.1367"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.911,
        "density": "1.0094"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.884,
        "density": "1.0168"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.864,
        "density": "1.0244"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.847,
        "density": "1.0319"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.832,
        "density": "1.0434"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.819,
        "density": "1.0511"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.808,
        "density": "1.0589"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.798,
        "density": "1.0668"
    },
    {
        "salt": "NaH2PO4",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.789,
        "density": "1.0747"
    },
    {
        "salt": "NaH2PO4",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.78,
        "density": "1.0826"
    },
    {
        "salt": "NaH2PO4",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.765,
        "density": "1.0988"
    },
    {
        "salt": "NaH2PO4",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.751,
        "density": "1.107"
    },
    {
        "salt": "NaH2PO4",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.739,
        "density": "1.1236"
    },
    {
        "salt": "NaH2PO4",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.729,
        "density": "1.1404"
    },
    {
        "salt": "NaH2PO4",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.721,
        "density": "1.149"
    },
    {
        "salt": "NaH2PO4",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.705,
        "density": "1.1931"
    },
    {
        "salt": "NaH2PO4",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.696,
        "density": "1.2113"
    },
    {
        "salt": "NaH2PO4",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.691,
        "density": "1.2488"
    },
    {
        "salt": "NaH2PO4",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 0.691,
        "density": "1.2879"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.793,
        "density": "1.00975"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.753,
        "density": "1.02593"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.725,
        "density": "1.03418"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.705,
        "density": "1.04617"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.69,
        "density": "1.05723"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.678,
        "density": "1.0713"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.667,
        "density": "1.0789"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.658,
        "density": "1.0905"
    },
    {
        "salt": "Na2SO4",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.65,
        "density": "1.10452"
    },
    {
        "salt": "Na2SO4",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.642,
        "density": "1.1122"
    },
    {
        "salt": "Na2SO4",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 0.631,
        "density": "1.12871"
    },
    {
        "salt": "Na2SO4",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 0.625,
        "density": "1.15369"
    },
    {
        "salt": "Na2SO4",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 0.621,
        "density": "1.17482"
    },
    {
        "salt": "Na2SO4",
        "molality": 1.8,
        "molarity": 1.802,
        "num_particles": 120,
        "osmotic_coefficient": 0.62,
        "density": "1.1907"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.767,
        "density": "1.0042"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.731,
        "density": "1.016"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.707,
        "density": "1.01968"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.69,
        "density": "1.0279"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.677,
        "density": "1.0338"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.667,
        "density": "1.042"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.658,
        "density": "1.0456"
    },
    {
        "salt": "NH42SO4",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.652,
        "density": "1.0574"
    },
    {
        "salt": "MgCl",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.861,
        "density": "1.0062"
    },
    {
        "salt": "MgCl",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.877,
        "density": "1.0144"
    },
    {
        "salt": "MgCl",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.895,
        "density": "1.0226"
    },
    {
        "salt": "MgCl",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.919,
        "density": "1.0309"
    },
    {
        "salt": "MgCl",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.947,
        "density": "1.0394"
    },
    {
        "salt": "MgCl",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.976,
        "density": "1.04629"
    },
    {
        "salt": "MgCl",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 1.004,
        "density": "1.0479"
    },
    {
        "salt": "MgCl",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 1.036,
        "density": "1.0564"
    },
    {
        "salt": "MgCl",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 1.071,
        "density": "1.0651"
    },
    {
        "salt": "MgCl",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 1.108,
        "density": "1.0738"
    },
    {
        "salt": "MgCl",
        "molality": 1.2,
        "molarity": 1.201,
        "num_particles": 80,
        "osmotic_coefficient": 1.184,
        "density": "1.0826"
    },
    {
        "salt": "MgCl",
        "molality": 1.4,
        "molarity": 1.401,
        "num_particles": 94,
        "osmotic_coefficient": 1.264,
        "density": "1.1005"
    },
    {
        "salt": "MgCl",
        "molality": 1.6,
        "molarity": 1.601,
        "num_particles": 107,
        "osmotic_coefficient": 1.347,
        "density": "1.11047"
    },
    {
        "salt": "MgCl",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 1.523,
        "density": "1.1372"
    },
    {
        "salt": "MgCl",
        "molality": 2.5,
        "molarity": 2.5,
        "num_particles": 167,
        "osmotic_coefficient": 1.762,
        "density": "1.240"
    },
    {
        "salt": "MgCl",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 2.01,
        "density": "1.1938"
    },
    {
        "salt": "MgCl",
        "molality": 3.5,
        "molarity": 4.003,
        "num_particles": 233,
        "osmotic_coefficient": 2.264,
        "density": "1.333"
    },
    {
        "salt": "MgCl",
        "molality": 4.0,
        "molarity": 4.003,
        "num_particles": 267,
        "osmotic_coefficient": 2.521,
        "density": "1.2555"
    },
    {
        "salt": "MgCl",
        "molality": 4.5,
        "molarity": 4.504,
        "num_particles": 300,
        "osmotic_coefficient": 2.783,
        "density": "1.2763"
    },
    {
        "salt": "NH42HP",
        "molality": 0.1,
        "molarity": 0.1,
        "num_particles": 7,
        "osmotic_coefficient": 0.89,
        "density": "1.004"
    },
    {
        "salt": "NH42HP",
        "molality": 0.2,
        "molarity": 0.2,
        "num_particles": 14,
        "osmotic_coefficient": 0.847,
        "density": "1.01"
    },
    {
        "salt": "NH42HP",
        "molality": 0.3,
        "molarity": 0.3,
        "num_particles": 20,
        "osmotic_coefficient": 0.816,
        "density": "1.017"
    },
    {
        "salt": "NH42HP",
        "molality": 0.4,
        "molarity": 0.4,
        "num_particles": 27,
        "osmotic_coefficient": 0.79,
        "density": "1.023"
    },
    {
        "salt": "NH42HP",
        "molality": 0.5,
        "molarity": 0.5,
        "num_particles": 34,
        "osmotic_coefficient": 0.77,
        "density": "1.029"
    },
    {
        "salt": "NH42HP",
        "molality": 0.6,
        "molarity": 0.601,
        "num_particles": 40,
        "osmotic_coefficient": 0.748,
        "density": "1.035"
    },
    {
        "salt": "NH42HP",
        "molality": 0.7,
        "molarity": 0.701,
        "num_particles": 47,
        "osmotic_coefficient": 0.738,
        "density": "1.041"
    },
    {
        "salt": "NH42HP",
        "molality": 0.8,
        "molarity": 0.801,
        "num_particles": 54,
        "osmotic_coefficient": 0.72,
        "density": "1.047"
    },
    {
        "salt": "NH42HP",
        "molality": 0.9,
        "molarity": 0.901,
        "num_particles": 60,
        "osmotic_coefficient": 0.71,
        "density": "1.05"
    },
    {
        "salt": "NH42HP",
        "molality": 1.0,
        "molarity": 1.001,
        "num_particles": 67,
        "osmotic_coefficient": 0.701,
        "density": "1.056"
    },
    {
        "salt": "NH42HP",
        "molality": 1.5,
        "molarity": 1.501,
        "num_particles": 100,
        "osmotic_coefficient": 0.658,
        "density": "1.081"
    },
    {
        "salt": "NH42HP",
        "molality": 2.0,
        "molarity": 2.002,
        "num_particles": 134,
        "osmotic_coefficient": 0.63,
        "density": "1.106"
    },
    {
        "salt": "NH42HP",
        "molality": 2.5,
        "molarity": 2.502,
        "num_particles": 167,
        "osmotic_coefficient": 0.608,
        "density": "1.126"
    },
    {
        "salt": "NH42HP",
        "molality": 3.0,
        "molarity": 3.003,
        "num_particles": 200,
        "osmotic_coefficient": 0.587,
        "density": "1.144"
    },
    {
        "salt": "NH42HP",
        "molality": 3.5,
        "molarity": 3.503,
        "num_particles": 234,
        "osmotic_coefficient": 0.573,
        "density": "1.159"
    }
]
