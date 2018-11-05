import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import flask
import json
from flask_cors import CORS
import os

app = dash.Dash('chembl-explorer')
server = app.server

df = pd.read_csv('lle_data.csv').drop(['target_chemblid.1', 'molregno.1'],
                                      axis=1)

axes_opts = [{'label': i, 'value': i} for i in df.columns if i not in ['target_chemblid','cmpd_chemblid',
                                                                       'full_molformula', 'target','molecular_species']]

axes_opts.append({'label':'rank X * Y','value':'rank'})

if 'DYNO' in os.environ:
    app.scripts.append_script({
        'external_url': 'https://cdn.rawgit.com/chriddyp/ca0d8f02a1659981a0ea7f013a378bbd/raw/e79f3f789517deec58f41251f7dbb6bee72c44ab/plotly_ga.js'
    })


def add_markers(figure_data, molecules, plot_type='scatter3d'):
    indices = []
    drug_data = figure_data[0]
    for m in molecules:
        hover_text = drug_data['text']
        for i in range(len(hover_text)):
            if m == hover_text[i]:
                indices.append(i)

    if plot_type == 'histogram2d':
        plot_type = 'scatter'

    traces = []
    for point_number in indices:
        trace = dict(
            x=[drug_data['x'][point_number]],
            y=[drug_data['y'][point_number]],
            marker=dict(
                color='red',
                size=16,
                opacity=0.6,
                symbol='cross'
            ),
            type=plot_type
        )

        if plot_type == 'scatter3d':
            trace['z'] = [drug_data['z'][point_number]]

        traces.append(trace)

    return traces


BACKGROUND = 'rgb(230, 230, 230)'

COLORSCALE = [[0, "rgb(244,236,21)"], [0.3, "rgb(249,210,41)"], [0.4, "rgb(134,191,118)"],
              [0.5, "rgb(37,180,167)"], [0.65, "rgb(17,123,215)"], [1, "rgb(54,50,153)"]]


def scatter_plot_3d(
        target_df=None,
        target=['Dihydrofolate reductase'],
        x_type='le',
        y_type='lle',
        z_type='rank',
        size_type='max_phase',
        color_type='rank',
        plot_type='scatter3d',
        markers=[]):
    def axis_template_3d(title, type='linear'):
        return dict(
            showbackground=True,
            backgroundcolor=BACKGROUND,
            gridcolor='rgb(255, 255, 255)',
            title=title,
            type=type,
            zerolinecolor='rgb(255, 255, 255)'
        )

    def axis_template_2d(title):
        return dict(
            xgap=10, ygap=10,
            backgroundcolor=BACKGROUND,
            gridcolor='rgb(255, 255, 255)',
            title=title,
            zerolinecolor='rgb(255, 255, 255)',
            color='#444'
        )

    def blackout_axis(axis):
        axis['showgrid'] = False
        axis['zeroline'] = False
        axis['color'] = 'white'
        return axis

    if target_df is None:
        target_df = df[df['target'] == target[0]]

    if z_type == 'rank':
        target_df['combined'] =  target_df[y_type] * target_df[x_type]
        target_df = target_df.sort_values(by=['combined'])
        z = [a for a in xrange(0, len(target_df['combined'].tolist()))]

    else:
        z = [a for a in target_df[z_type].tolist()]

    x = [b for b in target_df[x_type].tolist()]
    y = [c for c in target_df[y_type].tolist()]

    size = [(n + 1) * 200 for n in target_df[size_type].tolist()]
    # drugs = target_df[target_df['max_phase'] ==4]
    #
    # markers = drugs['cmpd_chemblid'].tolist()


    # if color_type == 'rank':
    #     try:
    #         color = [a for a in xrange(0, len(target_df['combined'].tolist()))]
    #     except KeyError:
    #         target_df['combined'] = target_df[y_type] * target_df[x_type]
    #         color = [a for a in xrange(0, len(target_df['combined'].tolist()))]
    # else:
    #     color = [d for d in target_df[color_type].tolist()]
    color = z
    xlabel = x_type
    ylabel = y_type
    zlabel = z_type

    data = [dict(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            colorscale='Viridis',
            # colorbar=dict(title="Molecular<br>Weight"),
            showscale=False,
            line=dict(color='#444'),
            sizeref=45,
            sizemode='diameter',
            opacity=0.7,
            size=size,
            color=color,
        ),
        text=target_df['cmpd_chemblid'].tolist(),
        type=plot_type,
    )]

    layout = dict(
        font=dict(family='Raleway'),
        hovermode='closest',
        margin=dict(r=20, t=0, l=0, b=0),
        showlegend=False,

        scene=dict(
            xaxis=axis_template_3d(xlabel),
            yaxis=axis_template_3d(ylabel),
            zaxis=axis_template_3d(zlabel),
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0.08, y=2.2, z=0.08)
            )
        )
    )

    if plot_type in ['histogram2d', 'scatter']:
        layout['xaxis'] = axis_template_2d(xlabel)
        layout['yaxis'] = axis_template_2d(ylabel)
        layout['plot_bgcolor'] = BACKGROUND
        layout['paper_bgcolor'] = BACKGROUND
        del layout['scene']
        del data[0]['z']

    if plot_type == 'histogram2d':
        # Scatter plot overlay on 2d Histogram
        data[0]['type'] = 'scatter'
        data.append(dict(
            x=x,
            y=y,
            type='histogram2d',
            colorscale='Viridis',#'Greys',
            showscale=False
        ))
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    if len(markers) > 0:
        data = data + add_markers(data, markers, plot_type=plot_type)

    return dict(data=data, layout=layout)


FIGURE = scatter_plot_3d()
DRUG_IMG = "https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL119"

app.layout = html.Div([
    # Row 1: Header and Intro text

    html.Div([
        html.H1('ChEMBL Explorer')
    ], className='row twelve columns', style={'position': 'relative', 'right': '15px'}),

    html.Div([
        html.Div([
            html.Div([
                html.P('SELECT a target in the dropdown to display its compounds in the graph'),
                html.P('CLICK on a compound in the graph to see its structure to the left.'),
                html.P('SELECT new axes from the dropdowns below. In all plots, the z axis will determine the colour'),
                html.P('Phase 4 compounds are marked with a red cross')

            ], style={'margin-left': '10px'}),
            dcc.Dropdown(id='chem_dropdown',
                         multi=True,
                         value=['Dihydrofolate reductase'],
                         options=[{'label': i, 'value': i} for i in df['target'].unique()]),
        ], className='twelve columns')

    ], className='row'),

    # Row 2: Hover Panel and Graph

    html.Div([
        html.Div([

            # html.Img(id='chem_img', src=DRUG_IMG, style=dict(width='100%')),
            html.Iframe(
                # enable all sandbox features
                # see https://developer.mozilla.org/en-US/docs/Web/HTML/Element/iframe
                # this prevents javascript from running inside the iframe
                # and other things security reasons
                id='chem_img',
                src="https://www.ebi.ac.uk/chembl/beta/embed/#mini_report_card/Compound/CHEMBL34259",
                style={'width': '100%', 'border': 'none'}
            ),

            html.Br(),




            dcc.Dropdown(id='x_dropdown',
                         multi=False,
                         value='lle',
                         options=axes_opts ),
            html.P("Select X axis"),



            dcc.Dropdown(id='y_dropdown',
                         multi=False,
                         value='le',
                         options=axes_opts ),
            html.P("Select Y axis"),

            dcc.Dropdown(id='z_dropdown',
                         multi=False,
                         value='rank',
                         options=axes_opts ),
            html.P("Select Z axis"),

        ], className='four columns'),

        html.Div([

            dcc.RadioItems(
                id='charts_radio',
                options=[
                    dict(label='3D Scatter', value='scatter3d'),
                    dict(label='2D Scatter', value='scatter'),
                    dict(label='2D Histogram', value='histogram2d'),
                ],
                labelStyle=dict(display='inline'),
                value='scatter3d'
            ),

            dcc.Graph(id='clickable-graph',
                      style=dict(width='700px', height='70vh'),
                      hoverData=dict(points=[dict(pointNumber=0)]),
                      figure=FIGURE),

        ], className='eight columns', style=dict(textAlign='center')),

    ], className='row'),

    # html.Div(id='signal', style={'display': 'none'}),

], className='container')


# # Store JSON in hidden div
# @app.callback(Output('signal', 'children'), [Input('chem_dropdown', 'value')])
# def global_store(n, value):
#     print('Computing value with {}'.format(value))
#     target_df = df[df['target'] ==value]
#     return target_df.to_json()

@app.callback(
    Output('clickable-graph', 'figure'),
    [Input('chem_dropdown', 'value'),
     Input('charts_radio', 'value'),
     Input('x_dropdown', 'value'),
     Input('y_dropdown', 'value'),
     Input('z_dropdown', 'value'),
     ])
def highlight_molecule(chem_dropdown_values, plot_type, x_type, y_type, z_type):
    return scatter_plot_3d(target=chem_dropdown_values, plot_type=plot_type, x_type=x_type, y_type=y_type, z_type=z_type)


# @app.callback(
#     State('signal','value'))
# def dfRowFromHover(hoverData, value):
#     ''' Returns row for hover point as a Pandas Series '''
#     print hoverData
#     print value
#     if hoverData is not None:
#         if 'points' in hoverData:
#             firstPoint = hoverData['points'][0]
#             if 'pointNumber' in firstPoint:
#                 point_number = firstPoint['pointNumber']
#                 molecule_name = str(FIGURE['data'][0]['text'][point_number]).strip()
#                 return df.loc[df['cmpd_chemblid'] == molecule_name]
#     return pd.Series()

#
# @app.callback(
#     Output('chem_name', 'children'),
#     [Input('clickable-graph', 'hoverData')])
# def return_molecule_name(hoverData):
#     if hoverData is not None:
#         if 'points' in hoverData:
#             firstPoint = hoverData['points'][0]
#             if 'pointNumber' in firstPoint:
#                 point_number = firstPoint['pointNumber']
#                 molecule_name = str(FIGURE['data'][0]['text'][point_number]).strip()
#                 return molecule_name




@app.callback(
    Output('chem_img', 'src'),
    [Input('clickable-graph', 'clickData')])
def display_image(hoverData):
    cmpd_id = hoverData['points'][0]['text']

    img_src = 'https://www.ebi.ac.uk/chembl/beta/embed/#mini_report_card/Compound/{}'.format(cmpd_id)
    return img_src


external_css = ["https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                "//fonts.googleapis.com/css?family=Raleway:400,300,600",
                "//fonts.googleapis.com/css?family=Dosis:Medium",
                "https://cdn.rawgit.com/plotly/dash-app-stylesheets/0e463810ed36927caf20372b6411690692f94819/dash-drug-discovery-demo-stylesheet.css"]

for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server()
