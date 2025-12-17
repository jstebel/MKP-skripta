# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "marimo",
#     "matplotlib==3.10.6",
#     "numpy==2.2.6",
# ]
# ///

import marimo

__generated_with = "0.18.0"
app = marimo.App(
    width="medium",
    css_file="/usr/local/_marimo/custom.css",
    auto_download=["html"],
)


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import numpy as np
    return np, plt


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # Plotting
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    marimo supports several popular plotting libraries, including matplotlib,
    plotly, seaborn, and altair.

    This tutorial gives examples using matplotlib; other libraries are
    used similarly.
    """)
    return


@app.cell
def _(np):
    xyz =1
    print(xyz)
    a = np.arange(10)
    print(a)
    return (xyz,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Matplotlib
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    To show a plot, include it in the last expression of a cell (just
    like any other output).

    ```python3
    # create the plot in the last line of the cell
    import matplotlib.pyplot as plt
    plt.plot([1, 2])
    ```
    """)
    return


@app.cell
def _(plt):
    plt.plot([1, 2])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ```python3
    # create a plot
    plt.plot([1, 2])
    # ... do some work ...
    # make plt.gca() the last line of the cell
    plt.gca()
    ```
    """)
    return


@app.cell
def _(plt):
    plt.plot([1, 2])
    # ... do some work ...
    # make plt.gca() the last line of the cell
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo, plt_show_explainer):
    mo.accordion(plt_show_explainer)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    **A new figure every cell.** Every cell starts with an empty figure for
    the imperative `pyplot` API.

    ## Try some formulas

    $x^2$
    """)
    return


@app.cell
def _(np):
    x = np.linspace(start=-4, stop=4, num=100, dtype=float)
    return (x,)


@app.cell
def _(plt, x, xyz):
    plt.plot(x, x)
    plt.plot(x, x**2)
    plt.gca()
    print(xyz)
    return


@app.cell
def _(plt, x):
    plt.plot(x, x**3)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    To build a figure over multiple cells, use the object-oriented API and
    create your own axis:
    """)
    return


@app.cell
def _(plt, x):
    _, axis = plt.subplots()
    axis.plot(x, x)
    axis.plot(x, x**2)
    axis
    return (axis,)


@app.cell
def _(axis, x):
    axis.plot(x, x**3)
    axis
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Draw plots interactively

    Draw plots interactively by parametrizing them with UI elements.
    """)
    return


@app.cell
def _(mo):
    exponent = mo.ui.slider(1, 5, value=1, step=1, label='exponent')

    mo.md(
        f"""
        **Visualizing powers.**

        {exponent}
        """
    )
    return (exponent,)


@app.cell
def _(plt, x):
    def plot_power(exponent):
        plt.plot(x, x**exponent)
        return plt.gca()
    return (plot_power,)


@app.cell
def _(exponent, mo, plot_power):
    _tex = (
        f"$$f(x) = x^{exponent.value}$$" if exponent.value > 1 else "$$f(x) = x$$"
    )

    mo.md(
        f"""

        {_tex}

        {mo.as_html(plot_power(exponent.value))}
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Other libraries
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    marimo also supports these other plotting libraries:

    - Plotly
    - Seaborn
    - Altair

    Just output their figure objects as the last expression of a cell,
    or embed them in markdown with `mo.as_html`.

    If you would like another library to be integrated into marimo, please
    get in touch.
    """)
    return


@app.cell
def _():
    plt_show_explainer = {
        "Using `plt.show()`": """
        You can use `plt.show()` or `figure.show()` to display
        plots in the console area of a cell. Keep in mind that console
        outputs are not shown in the app view.
        """
    }
    return (plt_show_explainer,)


if __name__ == "__main__":
    app.run()
