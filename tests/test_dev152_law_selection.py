from pbuf.foundation.native_neighbor_law_discriminator import decide

def test_equivalence_class_is_not_called_unique():
    rows=[{"viable":True}]
    classes=[{"members":["S-C10","S-C12"],"invertible_mapping_proven":True}]
    result=decide(rows,classes)
    assert result["outcome"]=="PBUF_NATIVE_NEIGHBOR_CONSTITUTIVE_EQUIVALENCE_CLASS_ESTABLISHED"
    assert not result["unique"]
