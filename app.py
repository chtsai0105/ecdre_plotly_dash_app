import pandas as pd
from dash import Dash, Input, Output, callback, dcc, html
from plotly import graph_objects as go

amplicon_types = ("16S", "ITS")
phylo_ranks = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

metadata_df = pd.read_csv("metadata.csv", index_col=0)
taxa_dfs = {
    amplicon_type: pd.read_csv(file, index_col=0) for amplicon_type, file in zip(amplicon_types, ("taxa_16s.csv", "taxa_its.csv"))
}
otu_dfs = {
    amplicon_type: pd.read_csv(file, index_col=0) for amplicon_type, file in zip(amplicon_types, ("otu_16s.csv", "otu_its.csv"))
}
treatment_dict_CA = {
    "CCC": "Control",
    "MCC": "Mulch",
    "CGC": "Glyphosate",
    "CCH": "Humic acid",
    "MGC": "Mulch + Glyphosate",
    "MCH": "Mulch + Humic acid",
    "CGH": "Glyphosate + Humic acid",
    "MGH": "Mulch + Glyphosate + Humic acid",
}
treatment_dict_app = {v: k for k, v in treatment_dict_CA.items()}

problematic_data = ["CCCT1R1", "L-T1-CCH1", "L-T1-CCH4", "L-T1-MCC4", "R-T1-MGC1"]


app = Dash(__name__)
server = app.server

app.layout = [
    html.Div(children="ECDRE year 0 and year 1 comparison"),
    html.Hr(),
    html.Label("Amplicon type:"),
    dcc.RadioItems(options=amplicon_types, value="16S", id="amplicon-type"),
    html.Div(
        children=[
            html.Div(
                children=[
                    html.Label("Tissue type:"),
                    dcc.Dropdown(options=metadata_df["Tissue_type"].unique(), value="Bulk soil", id="tissue", clearable=False),
                ],
                style={"padding": 10, "flex": 1},
            ),
            html.Div(
                children=[
                    html.Label("Treatment:"),
                    dcc.Dropdown(options=list(treatment_dict_app.keys()), value="Control", id="treatment", clearable=False),
                ],
                style={"padding": 10, "flex": 1},
            ),
            html.Div(
                children=[
                    html.Label("Phylogeny rank:"),
                    dcc.Dropdown(options=phylo_ranks[1:], value="Phylum", id="rank", clearable=False),
                ],
                style={"padding": 10, "flex": 1},
            ),
        ],
        style={"display": "flex", "flexDirection": "row"},
    ),
    html.Div(children="Suspected problematic data are prefixed with *"),
    dcc.Graph(figure={}, id="percentage-stacked-bar"),
]


def get_dataframe(amplicon_type, tissue, treatment, rank):
    otu_df = otu_dfs[amplicon_type]
    taxa_df = taxa_dfs[amplicon_type]
    metadata = metadata_df[metadata_df.index.isin(otu_df.columns)]
    metadata = metadata[(metadata["Treatment"] == treatment) & (metadata["Tissue_type"] == tissue)]
    otu_df = otu_df.loc[:, metadata.index]
    perc_df = otu_df / otu_df.sum() * 100
    perc_df = pd.merge(perc_df, taxa_df, left_index=True, right_index=True).groupby(by=list(phylo_ranks[: phylo_ranks.index(rank) + 1])).sum(numeric_only=True)
    perc_df = perc_df.reindex(perc_df.mean(axis=1).sort_values(ascending=False).index)
    perc_df = perc_df[perc_df.sum(axis=1) != 0]
    for sample in metadata[metadata.index.isin(problematic_data)].index:
        perc_df.rename({sample: f"*{sample}"}, axis=1, inplace=True)
        metadata.rename({sample: f"*{sample}"}, axis=0, inplace=True)
    return perc_df, metadata


# Add controls to build the interaction
@callback(
    Output("percentage-stacked-bar", "figure"),
    Input("amplicon-type", "value"),
    Input("tissue", "value"),
    Input("treatment", "value"),
    Input("rank", "value"),
)
def update_graph(amplicon_type, tissue, treatment, rank):
    perc_df, metadata = get_dataframe(amplicon_type, tissue, treatment_dict_app[treatment], rank)
    fig = go.Figure()
    for ids in perc_df.index:
        lineages = [taxon for taxon in ids if taxon != "Unclassified"]
        lineages = ["Unclassified"] if not lineages else lineages
        fig.add_trace(
            go.Bar(
                x=[[f"{x}, {y}" 
                    for x, y in zip(
                        metadata.loc[perc_df.columns, "Tissue_type"],
                        metadata.loc[perc_df.columns, "Year"])
                        ],
                    perc_df.columns],
                y=perc_df.loc[ids],
                name=lineages[-1],
                hovertemplate='<br>'.join(
                    [
                        '%{y:.4f}',
                        f'{";".join(lineages)}'
                    ]
                ),
            )
        )
    fig.update_layout(barmode="stack")
    return fig


if __name__ == "__main__":
    # app.run()
    app.run_server()
