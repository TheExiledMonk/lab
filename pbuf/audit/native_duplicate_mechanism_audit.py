def duplicates():
    return [
      {"concept":"trajectory","implementations":["excitation centroid history","pbuf.wl.trajectory_state geometric path","weak-lensing ray histories"],"classification":"DISTINCT_PHYSICS"},
      {"concept":"propagation","implementations":["np.roll excitation permutation","frame-aware link permutation","weak-lensing geometric propagation"],"classification":"SAME_PHYSICS_DIFFERENT_REPRESENTATION for first two; DISTINCT_PHYSICS for WL"},
      {"concept":"loading","implementations":["A8 c_state","Dev151 L profile"],"classification":"IMPLEMENTED_UNDER_DIFFERENT_NAME; no code derivation between them"},
      {"concept":"frame rotation","implementations":["frame_overlap","F02-F06 transport_map"],"classification":"SAME_PHYSICS_SAME_IMPLEMENTATION"},
      {"concept":"wavelength","implementations":["excitation FFT estimator","wl.native_wave_state metadata"],"classification":"DISTINCT_PHYSICS"},
      {"concept":"curvature","implementations":["WL geometric path curvature","claimed excitation trajectory curvature"],"classification":"HISTORICAL_ONLY for excitation; geometric implementation exists elsewhere"}]
