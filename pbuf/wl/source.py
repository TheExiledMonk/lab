"""Canonical benchmark-assisted source loading."""

import numpy as np

from pbuf.labs.foundation import current_native_five_cluster_observable_benchmark001 as CUR


def load_cluster_source(cluster: dict) -> dict:
    data = CUR.local_cluster(cluster)
    return {"cluster": cluster, "data": data, "rho3": np.asarray(data["rho3"])}
