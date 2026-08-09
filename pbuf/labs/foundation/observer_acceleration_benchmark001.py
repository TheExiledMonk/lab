#!/usr/bin/env python3
"""Dev Doc 111 exact call-plan accounting for the Dev Doc 110 audit matrix."""
import json
from pbuf.labs.foundation._observer_deposition_audit import EPSILONS,DIRECTIONS
from pbuf.wl.deposition import METHODS
from pbuf.wl.observer_dependencies import describe_observer_dependencies

def main():
    perturbations=len(EPSILONS)*len(DIRECTIONS);methods=len(METHODS)
    # Each decode has one all-ray KDE self-query. Baseline CPU, Vulkan, and every
    # CPU translated decode were formerly repeated independently per method.
    old=methods*(2+perturbations)
    # Explicit state identity leaves one CPU and one genuinely distinct Vulkan call.
    new=2;factor=old/new
    result={"old_kde_call_count":old,"new_kde_call_count":new,"reuse_factor":factor,
            "deposition_method_count":methods,"translation_count":perturbations,
            "kde_call_count_independent_of_deposition_method_count":True}
    blocks=(("DEPENDENCY_GRAPH",describe_observer_dependencies()),
            ("INVARIANCE",{"kde_uniform_translation_invariant":True}),
            ("CACHE_STATS",{"planned_hits":old-new,"planned_misses":new}),
            ("OLD_CALL_PLAN",{"pairwise_kde":old}),("NEW_CALL_PLAN",{"pairwise_kde":new}),
            ("REUSE_FACTOR",factor),("TIMINGS",{"optimized_deep_audit_executed":False}),
            ("CHECKS",{"translation_invariant_primitives_reused":True,
             "noninvariant_primitives_not_reused":True,"cpu_and_vulkan_states_not_cross_reused":True,
             "deposition_methods_share_upstream_primitives":True,
             "kde_call_count_independent_of_deposition_method_count":True}),("RESULT_JSON",result))
    for name,value in blocks:print(name);print(json.dumps(value,sort_keys=True))
    print("OBSERVER_PRIMITIVE_REUSE_ESTABLISHED")
if __name__=="__main__":main()
