use numpy::{PyReadonlyArray1, IntoPyArray};
use pyo3::prelude::*;
use quickseg_core::*;

#[pyfunction(signature = (
    values,
    penalty,
    segment_values=None,
    positions=None
))]
fn segment(
    py: Python,
    values: PyReadonlyArray1<f64>,
    penalty: f64,
    segment_values: Option<PyReadonlyArray1<f64>>,
    positions: Option<PyReadonlyArray1<f64>>,
) -> PyResult<Py<PyAny>> {

    let values = values.as_slice()?;
    assert!(!values.is_empty(), "Values must not be empty");

    let seg_values: Vec<f64> = match segment_values {
        Some(x) => x.as_slice()?.to_vec(),
        None => {
            let mut v = values.to_vec();
            v.sort_by(|a, b| a.total_cmp(b));
            v.dedup();
            v
        }
    };

    let (out_index, out_values) =
        segment_chromosome(values, &seg_values, penalty);

    let seg_end: Vec<usize> = out_index
        .iter()
        .skip(1)
        .map(|&x| x - 1)
        .chain(std::iter::once(values.len() - 1))
        .collect();

    if let Some(pos) = positions {
        let pos = pos.as_slice()?;

        let seg_begin: Vec<f64> =
            out_index.iter().map(|&i| pos[i]).collect();

        let seg_end_pos: Vec<f64> =
            seg_end.iter().map(|&i| pos[i]).collect();

        Ok((
            seg_begin.into_pyarray(py),
            seg_end_pos.into_pyarray(py),
            out_values.into_pyarray(py),
        )
            .into_pyobject(py)?
            .into())
    } else {
        Ok((
            out_index.into_pyarray(py),
            seg_end.into_pyarray(py),
            out_values.into_pyarray(py),
        )
            .into_pyobject(py)?
            .into())
    }
}

#[pymodule]
fn quickseg(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(segment, m)?)?;
    Ok(())
}
