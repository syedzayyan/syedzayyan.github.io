+++
title="Temporal Point Processes and Screaming Neighbours"
date="2026-02-11T14:00:01.000Z"

[taxonomies] 
tags = ["tutorial?"]
+++

I have been investigating into time-series models and landed onto Temporal Point Process one fine day. This world mostly concerts itself with finance and I am not a finance person or a statistics person. This is a tutorial or more appropiately my thoughts on how all of this works. Before I move along here is the problem we will be working with today.

### Problem

The issue at hand is my neighbour screams loudly everyday. Once only in a day but almost all days of the year. To simulate my neighbour doing professional neighbouring I am using a **Hawkes Process**. I will be using **Recurrent Marked Temporal Processes** to predict when my neighbour is going to scream again next, so that I am prepare.

```python {.marimo}
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
import torch.nn.functional as F
```

```python {.marimo}
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")
```

### SIMULATE SCREAM DATA

Before I explain what Hawkes Processes are, some background information is in order.

**Stochastic processes** are defined as (often infinite) collections of random variables, indexed by some time or space. In our case it is time, formalizing a model for the “discrete-time” data above, we could write {Xt}∞t=0,t∈Z+. A random variate, or a realization of the process is the entire trajectory determined by values taken by all Xt (often part of which we observe). In other words the scream, they are random, hence random variables. ⭐⭐⭐ *statistics*. Here Xt make up the collection while Z+ is the index set. The specific dependence (or rather, independence) structure, and other parametric assumptions of relationships among {Xt} determine the stochastic process.

A **Counting processes** is a stochastic process { N ( t ) , t ≥ 0 } {\displaystyle \{N(t),t\geq 0\}} with values that are non-negative, integer, and non-decreasing. In other words we count the number of screams in a day.

A **Poisson Process** is a counting process with quirks. Quirks being:

- We define a Poisson process with a function λ(t)>0,∀t
- Say we have two intervals, A,B⊂R The number of occurrences in these intervals will be Poisson distributed with N(A)∼∫Aλ(t)dt, and N(B)∼∫Bλ(t)dt
- Most importantly, N(A),N(B) are independent variables for all A∩B=∅
- Higher intensity functions λ(t)as expected, are associated with higher probabilities of event occurrences.

Imagine the neighbour screams because they are arguing. Arguments rarely end in one shout. One scream makes another scream more likely shortly after. But if a few calm days pass, the probability resets to baseline.

**Self-Exciting Processes** Now we allow events to influence future events. One scream increases the probability of another scream shortly after.But that influence decays over time. That’s called self-excitation. If neighbour screams because they are arguing:

- First scream → increases probability of second scream.
- After a few calm days → back to normal baseline.

Memory exists, but it fades.

##### Finally the, Hawkes Processes

The (univariate) Hawkes process is defined by the conditional intensity function

λ∗(t)=μ+∑ti<tφ(t−ti).

At any moment t, the conditional intensity function is at least μ>0, the background intensity. However, it also depends linearly on effects of events that have occurred before time t. Namely, this dependence is through a triggering kernel function φ(.), a function of the delay t−ti between the current time and the timestamp of the previous event. Note that φ is nonnegative (φ(x)≥0,∀x≥0 and causal φ(x)=0,∀x<0. It is usually a monotonically decreasing function (such as exponential decay, or power-law decay).

Thinking the other way around, the function can be interpreted as follows. Everytime neighbour screams, they'll scream more in the next days. *Self Exciting* called bursts. Probably, if they are sad, they'll remain sad for a while?

```python {.marimo}
def simulate_hawkes(T=365.0, mu=0.4, alpha=0.8, beta=1.5, seed=42):
    np.random.seed(seed)
    t = 0.0
    events = []

    while t < T:
        if len(events) == 0:
            lambda_bar = mu
        else:
            lambda_bar = mu + alpha * np.sum(np.exp(-beta * (t - np.array(events))))

        u = np.random.rand()
        w = -np.log(u) / lambda_bar
        t += w

        if t >= T:
            break

        lambda_t = mu + alpha * np.sum(np.exp(-beta * (t - np.array(events))))
        d = np.random.rand()

        if d * lambda_bar <= lambda_t:
            events.append(t)

    return np.array(events)

# Generate full data
true_times = simulate_hawkes()
gap_start, gap_end = 100.0, 130.0  # 30-day logger gap
latent_times = true_times[(true_times >= gap_start) & (true_times < gap_end)]
latent_count = len(latent_times)
obs_times = np.concatenate([
    true_times[true_times < gap_start],
    true_times[true_times >= gap_end]
])

print(f"True screams: {len(true_times)}, Observed: {len(obs_times)}, Latent in gap: {latent_count}")
```

Yes my neighbour seems to scream a lot

```python {.marimo}
plt.figure(figsize=(12, 3))

# Event (raster) plot
plt.plot(true_times, np.zeros_like(true_times), '|', markersize=12)

plt.xlabel("Time (days)")
plt.yticks([])
plt.title("Simulated Neighbour Screams (Hawkes Process)")

plt.tight_layout()
plt.show()
```

### PARAMETRIC HAWKES BASELINE

```python {.marimo}
class Hawkes(nn.Module):
    def __init__(self):
        super().__init__()
        self.mu = nn.Parameter(torch.tensor(0.3))
        self.alpha = nn.Parameter(torch.tensor(0.5))
        self.beta = nn.Parameter(torch.tensor(1.0))

    def intensity(self, t, history):
        past = history[history < t]
        if len(past) == 0:
            return self.mu
        return self.mu + self.alpha * torch.sum(
            torch.exp(-self.beta * (t - past))
        )

    def integral(self, a, b, history):
        mu_part = self.mu * (b - a)
        past = history[history < b]

        if len(past) == 0:
            return mu_part

        decay_a = torch.sum(torch.exp(-self.beta * (a - past)))
        decay_b = torch.sum(torch.exp(-self.beta * (b - past)))
        exc = (self.alpha / self.beta) * (decay_a - decay_b)

        return mu_part + exc

    def log_likelihood(self, times):
        ll = 0.0
        prev_t = torch.tensor(0.0, device=times.device)

        for i, t in enumerate(times):
            lamb = self.intensity(t, times[:i])
            integ = self.integral(prev_t, t, times[:i])
            ll += torch.log(lamb + 1e-8) - integ
            prev_t = t

        # subtract integral to T
        T = times[-1]
        ll -= self.integral(prev_t, T, times)

        return ll


def train_hawkes(model, times, epochs=400):
    optimizer = optim.Adam(model.parameters(), lr=0.05)
    times_t = torch.tensor(times, dtype=torch.float32, device=device)

    for e in range(epochs):
        optimizer.zero_grad()
        ll = model.log_likelihood(times_t)
        loss = -ll
        loss.backward()
        optimizer.step()

        if e % 100 == 0:
            print("Hawkes NLL:", -ll.item())

    br = model.alpha.item() / model.beta.item()
    print("Branching ratio α/β:", br)
```

### RMTPP MODEL (Recurrent Marked TPP, simplified univariate)
<!---->
So far, Hawkes assumed memory decays exponentially. But what if your neighbour’s emotional state is not exponential? Instead of hard-coding the memory kernel, we can let a neural network learn it.

That is RMTPP.

We observe inter-event times:
$$
\tau_j = t_j - t_{j-1}
$$
These are fed into an RNN (typically an LSTM), producing hidden states:
$$
h_j = \mathrm{RNN}(h_{j-1}, \tau_j)
$$

The hidden state $h_j$ is a learned summary of the entire past. This replaces the explicit Hawkes kernel. The [original paper](https://www.kdd.org/kdd2016/papers/files/rpp1081-duA.pdf) is a way better resource. Please head over there if you want more insights.
<!---->
| Hawkes                    | RMTPP                        |
| ------------------------- | ---------------------------- |
| Exponential memory kernel | Learned memory               |
| Interpretable parameters  | Black box                    |
| Parametric                | Neural                       |
| Stable if α/β < 1         | Stability learned implicitly |

```python {.marimo}
class RMTPP(nn.Module):
    def __init__(self, hidden_dim=32):
        super().__init__()
        self.lstm = nn.LSTM(1, hidden_dim, batch_first=True)
        self.v = nn.Linear(hidden_dim, 1)
        self.w = nn.Parameter(torch.tensor(0.1))
        self.b = nn.Parameter(torch.tensor(0.0))

    def forward_history(self, deltas):
        deltas = deltas.unsqueeze(-1)
        out, _ = self.lstm(deltas)
        return out

    def intensity(self, h, tau):
        return torch.exp(self.v(h) + self.w * tau + self.b)

    def log_density(self, h, tau):
        a = self.v(h) + self.b
        w = self.w

        lamb = torch.exp(a + w * tau)
        integral = (torch.exp(a) / w) * (torch.exp(w * tau) - 1)

        return torch.log(lamb + 1e-8) - integral

    def log_likelihood(self, deltas):
        h_all = self.forward_history(deltas.unsqueeze(0))[0]

        ll = 0.0
        for j in range(1, len(deltas)):
            h_prev = h_all[j - 1]
            tau = deltas[j]
            ll += self.log_density(h_prev, tau)

        return ll
```

### TRAINING (Hawkes + RMTPP)

```python {.marimo}
def train_rmtpp(model, times, epochs=600):
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    deltas = np.diff(np.concatenate([[0.0], times]))
    deltas_t = torch.tensor(deltas, dtype=torch.float32, device=device)

    for e in range(epochs):
        optimizer.zero_grad()
        ll = model.log_likelihood(deltas_t)
        loss = -ll
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

        if e % 150 == 0:
            print("RMTPP NLL:", -ll.item())

def sample_rmtpp_next(model, times):
    deltas = np.diff(np.concatenate([[0.0], times]))
    deltas_t = torch.tensor(deltas, dtype=torch.float32, device=device)

    h_all = model.forward_history(deltas_t.unsqueeze(0))[0]
    h_last = h_all[-1]

    a = model.v(h_last) + model.b
    w = model.w

    u = torch.rand(1, device=device)

    tau = (1.0 / w) * torch.log(
        1 - (w / torch.exp(a)) * torch.log(1 - u)
    )

    return tau.item()
```

```python {.marimo}
hawkes = Hawkes().to(device)
train_hawkes(hawkes, true_times)

rmtpp = RMTPP().to(device)
train_rmtpp(rmtpp, true_times)

next_tau = sample_rmtpp_next(rmtpp, true_times)
print("Predicted next scream in ~", next_tau, "days")
```

### PLOT RESULTS

```python {.marimo}
t_eval = np.linspace(true_times[-1], true_times[-1] + 10, 200)
history_t = torch.tensor(true_times, dtype=torch.float32)

hawkes_int = [
    hawkes.intensity(torch.tensor(t), history_t).item()
    for t in t_eval
]

plt.figure(figsize=(10,5))
plt.plot(t_eval, hawkes_int)
plt.title("Hawkes Predicted Intensity")
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.show()
```

References:

- https://shchur.github.io/blog/2020/tpp1-conditional-intensity/
- https://hawkeslib.readthedocs.io/en/latest/tutorial.html#footnote-reference-1
- https://en.wikipedia.org/wiki/Counting_process
