#!/usr/bin/env python3
"""Structural verification of the hand-written ODE (henze2017.rhs).

Independent reference: evaluate the kinetic laws DIRECTLY from the raw
SBML XML (c2v MathML AST, parsed fresh here), combined with the reaction
stoichiometry read from the same XML. No intermediate JSON, no shared
code with henze2017.py. Agreement at 10,000 random nonnegative states
proves the hand-written rhs is 1:1 with the deposited network.
"""
import json, os, re, sys
import xml.etree.ElementTree as ET
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import henze2017 as H

SBML = 'http://www.sbml.org/sbml/level2/version4'
MATH = 'http://www.w3.org/1998/Math/MathML'
def q(t): return f'{{{SBML}}}{t}'
def mq(t): return f'{{{MATH}}}{t}'

tree = ET.parse(os.path.join(HERE, "Henze2017_original_SBML.xml"))
model = tree.getroot().find(q('model'))
spec = [s.get('id') for s in model.iter(q('species'))]
assert set(spec) == set(H.SPECIES), "species sets differ"
# xml order may differ from module order; build permutation xml -> module
perm = np.array([spec.index(s) for s in H.SPECIES])  # H.SPECIES[i] = xml spec[perm[i]]
params = {p.get('id'): float(p.get('value')) for p in model.iter(q('parameter'))}

# --- c2v AST evaluator (fresh implementation) -----------------------------
OPS = {'plus': lambda a, b: a + b, 'minus': lambda a, b: a - b,
       'times': lambda a, b: a * b, 'divide': lambda a, b: a / b,
       'power': lambda a, b: a ** b}

def eval_math(el, env):
    tag = el.tag
    if tag in (mq('ci'), mq('cn')):
        name = (el.text or '').strip()
        return float(name) if re.fullmatch(r'\d+(\.\d+)?', name) else env[name]
    if tag == mq('apply'):
        kids = list(el)
        op0 = kids[0]
        fname = op0.find(mq('name')).text if op0.tag == mq('apply') else op0.tag.split('}')[1]
        vals = [eval_math(a, env) for a in kids[1:]]
        if fname in OPS:
            acc = vals[0]
            for v in vals[1:]:
                acc = OPS[fname](acc, v)
            return acc
        raise ValueError(f"unsupported c2v function {fname}")
    if tag in (mq('math'), mq('semantics')):
        for c in el:
            r = eval_math(c, env)
            if r is not None:
                return r
        return None
    if tag == mq('annotation'):
        return None
    if tag == mq('annotation-xml'):
        for c in el:
            r = eval_math(c, env)
            if r is not None:
                return r
        return None
    raise ValueError(f"unsupported MathML element {tag}")

# --- reactions from raw XML ------------------------------------------------
reactions = []
for r in model.iter(q('reaction')):
    rec = {'reactants': {}, 'products': {}, 'law_el': None}
    for side, key in [(q('listOfReactants'), 'reactants'), (q('listOfProducts'), 'products')]:
        lst = r.find(side)
        if lst is None:
            continue
        for sp in lst.iter(q('speciesReference')):
            s = sp.get('species')
            st = int(sp.get('stoichiometry', '1'))
            rec[key][s] = rec[key].get(s, 0) + st
    kl = r.find(q('kineticLaw'))
    rec['law_el'] = kl.find(mq('math'))
    assert rec['law_el'] is not None, f"reaction {r.get('id')} has no math"
    reactions.append(rec)

# assignment rules (variable="..." attribute form)
assign_rules = {}
for r in model.iter(q('assignmentRule')):
    assign_rules[r.get('variable')] = r.find(mq('math'))

def env_at(y):
    env = {s: float(y[i]) for i, s in enumerate(spec)}
    env.update(params)
    env['compartment'] = 1.0
    for var, math_el in assign_rules.items():
        env[var] = eval_math(math_el, env)
    return env

def generic_rhs(y):
    env = env_at(y)
    d = np.zeros(len(spec))
    for rec in reactions:
        v = eval_math(rec['law_el'], env)
        for s, st in rec['reactants'].items():
            d[spec.index(s)] -= st * v
        for s, st in rec['products'].items():
            d[spec.index(s)] += st * v
    return d

# --- compare ----------------------------------------------------------------
rng = np.random.default_rng(42)
y0 = np.array([float(s.get('initialConcentration') or 0) for s in model.iter(q('species'))])
n_states = 10000
maxdiff, worst = 0.0, None
for i in range(n_states):
    scale = np.array([200.0 if s in ("KinU", "KinA") else 2.5 * max(y0[i], 0.05)
                      for i, s in enumerate(spec)])
    y = rng.random(len(spec)) * scale
    a = generic_rhs(y)                      # xml species order
    y_h = y[perm]                           # module species order
    b_h = H.rhs(0.0, y_h)
    b = np.empty(len(spec)); b[perm] = b_h  # back to xml order
    dd = float(np.max(np.abs(a - b)))
    if dd > maxdiff:
        maxdiff, worst = dd, i
print(f"states tested: {n_states}")
print(f"max |raw-SBML-generic - handwritten-rhs| = {maxdiff:.3e}")
if maxdiff < 1e-9:
    print("VERDICT: PASS - handwritten rhs is 1:1 with the deposited SBML network")
else:
    print("VERDICT: FAIL at state index", worst)
    sys.exit(1)
