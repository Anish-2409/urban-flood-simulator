"""Streamlit interactive dashboard for the 2D urban flood simulator.

Educational interface enabling users to:
1. Select topography (bowl, slope, random microtopography).
2. Configure storm events (constant vs. triangular design storm).
3. Compare baseline vs. grey drainage vs. green infrastructure scenarios.
4. Visualize water depth maps and mass balance conservation curves.
"""

import streamlit as st


def main() -> None:
    """Run Streamlit user interface."""
    st.set_page_config(page_title="Urban Flood Simulator", layout="wide")
    st.title("🌧️ Urban Flood Simulator (2D Waterlogging)")
    st.caption("Educational 2D numerical water balance model for engineering mathematics.")

    st.sidebar.header("Scenario Settings")
    scenario_choice = st.sidebar.selectbox(
        "Select Scenario",
        ["Baseline (Impervious Urban)", "Improved Drains", "Green Infrastructure (Sponge City)"],
    )

    st.sidebar.subheader("Rainfall Hyetograph")
    storm_type = st.sidebar.radio("Storm Profile", ["Constant", "Triangular"])
    intensity = st.sidebar.slider("Intensity / Peak (mm/hr)", min_value=10.0, max_value=120.0, value=50.0)
    duration = st.sidebar.slider("Duration (minutes)", min_value=10.0, max_value=120.0, value=30.0)

    st.info(
        f"Selected scenario: **{scenario_choice}** with {storm_type.lower()} storm "
        f"({intensity} mm/hr for {duration} min).\n\n"
        "Numerical simulation backend is currently configured as contract stubs. "
        "Run `pytest tests/` to verify contracts."
    )


if __name__ == "__main__":
    main()
