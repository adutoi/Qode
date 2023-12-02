#    (C) Copyright 2023 Anthony D. Dutoi
# 
#    This file is part of Qode.
# 
#    Qode is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    Qode is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with Qode.  If not, see <http://www.gnu.org/licenses/>.
#
import numpy

def one_electron(O_raw, p_blocks, q_blocks=None):
    norb_p, norb_q = O_raw.shape
    #
    if q_blocks is None:  q_blocks = p_blocks    # so by default, p_blocks and q_blocks are the same
    if p_blocks=="blocked":      p_blocks = [[p for p in range(norb_p)]]    # If p_blocks and q_blocks are not ...
    if p_blocks=="alternating":  p_blocks = [[p] for p in range(norb_p)]    # ... one of these special strings then
    if q_blocks=="blocked":      q_blocks = [[q for q in range(norb_q)]]    # ... they should be lists of lists containingg
    if q_blocks=="alternating":  q_blocks = [[q] for q in range(norb_q)]    # ... the indices of each spin block.
    if norb_p!=sum(len(block) for block in p_blocks):  raise ValueError
    if norb_q!=sum(len(block) for block in q_blocks):  raise ValueError
    #
    O = numpy.zeros((2*norb_p, 2*norb_q))
    #
    p0_dn = 0
    for p_block in p_blocks:
        p0_up = p0_dn + len(p_block)
        for p_,p in enumerate(p_block):
            q0_dn = 0
            for q_block in q_blocks:
                q0_up = q0_dn + len(q_block)
                for q_,q in enumerate(q_block):
                    O[p0_dn+p_, q0_dn+q_] = O_raw[p, q]
                    O[p0_up+p_, q0_up+q_] = O_raw[p, q]
                q0_dn += 2*len(q_block)
        p0_dn += 2*len(p_block)
    #
    return O

def two_electron(O_raw, p_blocks, q_blocks=None, r_blocks=None, s_blocks=None):
    norb_p, norb_q, norb_r, norb_s = O_raw.shape
    #
    if q_blocks is None:  q_blocks = p_blocks    # so by default, p_blocks and q_blocks are the same
    if r_blocks is None:  r_blocks = p_blocks    # so by default, p_blocks and r_blocks are the same
    if s_blocks is None:  s_blocks = q_blocks    # so by default, q_blocks and s_blocks are the same
    if p_blocks=="blocked":      p_blocks = [[p for p in range(norb_p)]]    # See comments in one_electron function ...
    if p_blocks=="alternating":  p_blocks = [[p] for p in range(norb_p)]
    if q_blocks=="blocked":      q_blocks = [[q for q in range(norb_q)]]
    if q_blocks=="alternating":  q_blocks = [[q] for q in range(norb_q)]
    if r_blocks=="blocked":      r_blocks = [[r for r in range(norb_r)]]
    if r_blocks=="alternating":  r_blocks = [[r] for r in range(norb_r)]
    if s_blocks=="blocked":      s_blocks = [[s for s in range(norb_s)]]
    if s_blocks=="alternating":  s_blocks = [[s] for s in range(norb_s)]
    if norb_p!=sum(len(block) for block in p_blocks):  raise ValueError
    if norb_q!=sum(len(block) for block in q_blocks):  raise ValueError
    if norb_r!=sum(len(block) for block in r_blocks):  raise ValueError
    if norb_s!=sum(len(block) for block in s_blocks):  raise ValueError
    #
    O = numpy.zeros((2*norb_p, 2*norb_q, 2*norb_r, 2*norb_s))
    #
    p0_dn = 0
    for p_block in p_blocks:
        p0_up = p0_dn + len(p_block)
        for p_,p in enumerate(p_block):
            q0_dn = 0
            for q_block in q_blocks:
                q0_up = q0_dn + len(q_block)
                for q_,q in enumerate(q_block):
                    r0_dn = 0
                    for r_block in r_blocks:
                        r0_up = r0_dn + len(r_block)
                        for r_,r in enumerate(r_block):
                            s0_dn = 0
                            for s_block in s_blocks:
                                s0_up = s0_dn + len(s_block)
                                for s_,s in enumerate(s_block):
                                    O[p0_dn+p_, q0_dn+q_, r0_dn+r_, s0_dn+s_] = O_raw[p, q, r, s]
                                    O[p0_up+p_, q0_dn+q_, r0_up+r_, s0_dn+s_] = O_raw[p, q, r, s]
                                    O[p0_dn+p_, q0_up+q_, r0_dn+r_, s0_up+s_] = O_raw[p, q, r, s]
                                    O[p0_up+p_, q0_up+q_, r0_up+r_, s0_up+s_] = O_raw[p, q, r, s]
                                s0_dn += 2*len(s_block)
                        r0_dn += 2*len(r_block)
                q0_dn += 2*len(q_block)
        p0_dn += 2*len(p_block)
    #
    return O

def one_electron_blocked(O_raw):
    #raise NotImplementedError("This function deprecated.  Replace with spatial_to_spin.one_electron(<raw>, \"blocked\")")
    return one_electron(O_raw, "blocked")

def one_electron_alternating(O_raw):
    #raise NotImplementedError("This function deprecated.  Replace with spatial_to_spin.one_electron(<raw>, \"alternating\")")
    return one_electron(O_raw, "alternating")

def two_electron_blocked(O_raw):
    #raise NotImplementedError("This function deprecated.  Replace with spatial_to_spin.two_electron(<raw>, \"blocked\")")
    return two_electron(O_raw, "blocked")

def two_electron_alternating(O_raw):
    #raise NotImplementedError("This function deprecated.  Replace with spatial_to_spin.two_electron(<raw>, \"alternating\")")
    return two_electron(O_raw, "alternating")
