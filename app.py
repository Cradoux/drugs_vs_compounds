import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import flask
import json
from flask_cors import CORS
import os
import numpy as np

app = dash.Dash('chembl-explorer')
server = app.server

df = pd.read_csv('lle_data_single_prot_human.csv')
axes_opts = [{'label': i, 'value': i} for i in df.columns if i not in ['target_chemblid', 'cmpd_chemblid',
                                                                       'target']]

axes_opts.append({'label': 'rank X * Y', 'value': 'rank'})

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
            showlegend=False,
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

        target='Integrin alpha-V/beta-3_PROTEIN COMPLEX_Homo sapiens_104292',
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

    data = []

    target_df = df[df['target_id'] == target]

    if z_type == 'rank':
        target_df['combined'] = target_df[y_type] * target_df[x_type]
        target_df = target_df.sort_values(by=['combined'])
        z = [a for a in xrange(0, len(target_df['combined'].tolist()))]

    else:
        z = [a for a in target_df[z_type].tolist()]

    x = [b for b in target_df[x_type].tolist()]
    y = [c for c in target_df[y_type].tolist()]

    size = [(n + 1) * 300 for n in target_df[size_type].tolist()]
    drugs = target_df[target_df['max_phase'] == 4]
    markers = []
    markers = markers + drugs['cmpd_chemblid'].tolist()

    color = z
    xlabel = x_type
    ylabel = y_type
    zlabel = z_type

    data.append(dict(
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
    ))

    if len(markers) > 0:
        data = data + add_markers(data, markers, plot_type=plot_type)

    layout = dict(
        font=dict(family='Raleway'),
        hovermode='closest',
        margin=dict(r=20, t=0, l=0, b=0),
        showlegend=False,
        legend=dict(orientation="h"),

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
            colorscale='Greys',
            showscale=False
        ))
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    return dict(data=data, layout=layout)


def summary_plot(df=df,
                 x_type='le',
                 z_type=None,
                 size_type='max_phase',
                 color_type='rank',
                 plot_type='scatter',
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

    data = []

    drugs = df[df['max_phase'] == 4]
    not_drugs = df[df['max_phase'] < 4]


    drugs = drugs.groupby('target_id', as_index=False).mean()
    not_drugs = not_drugs.groupby('target_id', as_index=False).mean()
    not_drugs = not_drugs[not_drugs['target_id'].isin(drugs['target_id'].unique())]
    text = drugs['target_id'].tolist()


    x = [b for b in not_drugs[x_type].tolist()]
    y = [c for c in drugs[x_type].tolist()]

    xlabel = 'Mean ' + x_type + ' Discovery compounds'
    ylabel = 'Mean ' + x_type + ' Drugs'
    zlabel = z_type





    data.append(dict(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            colorscale='Viridis',
            # colorbar=dict(title="Molecular<br>Weight"),
            showscale=False,
            line=dict(color='#444'),
            sizeref=45,
            sizemode='diameter',
            opacity=0.7,
        ),
        text=text,
        type=plot_type,
    ))

    markers = [markers]

    if len(markers) > 0:

        data = data + add_markers(data, markers, plot_type=plot_type)

    layout = dict(
        font=dict(family='Raleway'),
        hovermode='closest',
        margin=dict(r=20, t=0, l=40, b=40),
        showlegend=False,
        legend=dict(orientation="h"),

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

    if plot_type == 'histogram2d':
        # Scatter plot overlay on 2d Histogram
        data[0]['type'] = 'scatter'
        data.append(dict(
            x=x,
            y=y,
            type='histogram2d',
            colorscale='Greys',
            showscale=False
        ))
        layout['plot_bgcolor'] = 'black'
        layout['paper_bgcolor'] = 'black'
        layout['xaxis'] = blackout_axis(layout['xaxis'])
        layout['yaxis'] = blackout_axis(layout['yaxis'])
        layout['font']['color'] = 'white'

    max_val = max(max(x),max(y))

    line = np.linspace(0,max_val,5)

    data.append(dict(
        x = line,
        y = line,
        mode = 'line'
    ))

    return dict(data=data, layout=layout)


FIGURE = scatter_plot_3d()
SUMMARY = summary_plot()
DRUG_IMG = "https://www.ebi.ac.uk/chembl/api/data/image/CHEMBL119"

app.layout = html.Div([
    # Row 1: Header and Intro text

    html.Div([
        html.H1('ChEMBL Explorer', id='header')
    ], className='mui-row ', style={'position': 'relative'}),

    # html.Div([
    #     html.Div([
    #         # html.Div([
    #         #     html.P('SELECT one or more targets in the dropdown to display their compounds in the graph'),
    #         #     html.P('CLICK on a compound in the graph to see its structure to the left.'),
    #         #     html.P('SELECT new axes from the dropdowns below. In all plots, the z axis will determine the colour'),
    #         #     html.P('Phase 4 compounds are marked with a red cross')
    #         #
    #         # ], style={'margin-left': '10px'}),
    #         dcc.Dropdown(id='chem_dropdown',
    #                      multi=False,
    #                      value='Dihydrofolate reductase',
    #                      options=[{'label': i, 'value': i} for i in df['target'].unique()]),
    #     ], className='mui-col-md-12')
    #
    # ], className='mui-row'),

    # Row 2: Hover Panel and Graph

    html.Div([
        html.Div([

            html.Div([

                html.Label('Select Property One'),
                dcc.Dropdown(id='x_dropdown',
                             multi=False,
                             value='le',
                             options=axes_opts),
                dcc.Graph(id='target-x-graph',
                          style=dict(height='35vh'),
                          hoverData=dict(points=[dict(pointNumber=0)]),
                          figure=summary_plot()),
            ], className='mui-row'),
            html.Div([
                html.Label('Select Property Two'),
                dcc.Dropdown(id='y_dropdown',
                             multi=False,
                             value='lle',
                             options=axes_opts),
                dcc.Graph(id='target-y-graph',
                          style=dict(height='35vh'),
                          hoverData=dict(points=[dict(pointNumber=0)]),
                          figure=summary_plot()),

            ], className='mui-row'),

        ], className='mui-col-md-5 mui-col-xs-12 '),

        html.Div([

            html.Label('Select Property Three'),
            dcc.Dropdown(id='z_dropdown',
                         multi=False,
                         value='rank',
                         options=axes_opts),

            dcc.Graph(id='clickable-graph',
                      style=dict(height='50vh'),
                      hoverData=dict(points=[dict(pointNumber=0)]),
                      figure=FIGURE),

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

            html.Div([

            # html.Img(id='chem_img', src=DRUG_IMG, style=dict(width='100%')),
            html.Iframe(
                id='chem_img',
                src="https://www.ebi.ac.uk/chembl/beta/embed/#mini_report_card/Compound/CHEMBL1046",
                style={'width': '100%', 'height': '30vh', 'border': 'none'}
            ),

        ], style=dict(textAlign='center', height='100%')),

        ], className='mui-col-md-7 mui-col-xs-12 ', style=dict(textAlign='center')),



    ], className='mui-row ', style=dict(textAlign='center', height='80vh')),
    html.Div(id='signal', style={'display': 'none'}),

], className='mui-container mui-panel')


# Store target JSON in hidden div
@app.callback(Output('signal', 'children'),
              [#Input('chem_dropdown', 'value'),
               Input('charts_radio', 'value'),
               Input('x_dropdown', 'value'),
               Input('y_dropdown', 'value'),
               Input('z_dropdown', 'value'),
               Input('target-x-graph', 'clickData'),
               Input('target-y-graph', 'clickData')
               ])
def target_store( plot_type, x_type, y_type, z_type,target_x,target_y):
    try:
        t_y = target_y['points'][0]['text']
    except TypeError:
        t_y = None

    try:
        t_x = target_x['points'][0]['text']
    except TypeError:
        t_x = None

    j = {
        'target':'Plasminogen_SINGLE PROTEIN_Homo sapiens_12',
        'plot_type':plot_type,
        'x_type': x_type,
        'y_type':y_type,
        'z_type':z_type,
        'target_x':t_x,
        'target_y':t_y
    }

    if t_x is not None:
        j['target'] = j['target_x']
        j['target_y'] = j['target_x']

    if t_y is not None:
        j['target'] = j['target_y']
        j['target_x'] = j['target_y']


    return json.dumps(j)


@app.callback(
    Output('clickable-graph', 'figure'),
    [Input('signal', 'children')])
def update_compound_graph(j):
    j = json.loads(j)

    return scatter_plot_3d(target=j['target'], plot_type=j['plot_type'], x_type=j['x_type'], y_type=j['y_type'],
                           z_type=j['z_type'])


@app.callback(
    Output('target-x-graph', 'figure'),
    [Input('signal', 'children')])
def update_target_x(j):
    j = json.loads(j)

    return summary_plot(markers=j['target'],
                        x_type=j['x_type'],
                           z_type=j['z_type'])

@app.callback(
    Output('target-y-graph', 'figure'),
    [Input('signal', 'children')])
def update_target_y(j):
    j = json.loads(j)

    return summary_plot(markers=j['target'],
                        x_type=j['y_type'],
                           z_type=j['z_type'])

@app.callback(
    Output('header', 'children'),
    [Input('signal', 'children')])
def update_target_y(j):
    j = json.loads(j)

    return "ChEMBL Explorer - {}".format(j['target'])


@app.callback(
    Output('chem_img', 'src'),
    [Input('clickable-graph', 'clickData')])
def display_image(hoverData):
    cmpd_id = hoverData['points'][0]['text']
    img_src = 'https://www.ebi.ac.uk/chembl/beta/embed/#mini_report_card/Compound/{}'.format(cmpd_id)
    # img_src = 'https://www.ebi.ac.uk/chembl/beta/embed/#compound_report_card/{}/name_and_classification'.format(cmpd_id)

    return img_src


external_css = ["https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
                "//fonts.googleapis.com/css?family=Raleway:400,300,600",
                "//fonts.googleapis.com/css?family=Dosis:Medium",
                "//cdn.muicss.com/mui-0.9.38/js/mui.min.js",
                "//cdn.muicss.com/mui-0.9.38/css/mui.min.css",
                "https://codepen.io/mikesmith1611/pen/QOKgpG.css"
                ]
# "https://cdn.rawgit.com/plotly/dash-app-stylesheets/0e463810ed36927caf20372b6411690692f94819/dash-drug-discovery-demo-stylesheet.css"]

for css in external_css:
    app.css.append_css({"external_url": css})

if __name__ == '__main__':
    app.run_server()
