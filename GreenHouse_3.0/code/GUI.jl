using PlotlyJS
using Dash
using Dates
using JSON3

include("simulation_GUI.jl")

using .HarverstSim

# Simulated real-time data
global data_history = Dict(
    "time" => Float64[],
    "T1" => Float64[],
    "T2" => Float64[],
    "CO2" => Float64[],
    "VP" => Float64[],
    "A" => Float64[],
    "RH" => Float64[],
    "PAR" => Float64[],
    "Tp" => Float64[],
    "SW" => Float64[],
    "Tout" => Float64[],
    "Tsky" => Float64[],
    "Tsoil" => Float64[],
    "Wind" => Float64[],
    "RHout" => Float64[],
    "VPout" => Float64[]

)

function clean_dict(dict::Dict{String, Vector{Float64}})
    return Dict(k => filter(!isnan, v) for (k, v) in dict)
end

function updateData(dias, d)
    for key in keys(d)
        d[key] = round.(100*rand(144*dias);digits = 1)
    end 
    return d
end


function deep_collect(x)
    if x isa JSON3.Array
        return [deep_collect(e) for e in x]
    elseif x isa JSON3.Object
        return Dict(string(k) => deep_collect(v) for (k, v) in x)
    elseif x isa Dict
        return Dict(k => deep_collect(v) for (k, v) in x)
    else
        return x
    end
end

function create_text_block(value, label, color, left, top, arrow)
    html_div("$label: $(round(value, digits=1)) $arrow", style=Dict(
        :position => "absolute",
        :left => left,
        :top => top,
        :color => color,
        :fontFamily=>"Arial",
        #:fontWeight => "bold",
        :fontSize => "15px",
        :backgroundColor => "rgba(20,20,20)",
        :padding => "4px 8px",
        :borderRadius => "5px"
    ))
end

app = dash(external_stylesheets = ["assets/GUI.css"])
app.title = "Greenhouse simulator"
app.layout = html_div() do
    dcc_interval(id="interval", interval=200, n_intervals=0, disabled = true),
    html_h1("Greenhouse simulation", style = Dict("margin-top" => 50,:textAlign => "center")),
    html_div(id="container",
        children = 
            [  # Wrapper
            html_img(id="greenhouse-img", src="assets/greenhouse3.png"), 
            html_div(id = "overlay", children=[]),
               html_div(
            children = 
                [
                html_h3("Selecciona fecha de inicio"),    
                dcc_datepickersingle(
                    id="date-picker",
                    min_date_allowed = Date(1998, 7, 11),
                    max_date_allowed = Date(2018, 5, 1),
                    initial_visible_month= Date(2017, 8, 5),
                    date = Date(2017, 8, 25)
                ),
                html_h3("Días de cultivo"),
                dcc_input(id = "dias_sim",  value = 7, type = "number", min=1, max=120),
                html_br(),
                html_br(),
                html_button("Upload simulation", id="sim-button", n_clicks=0),
                dcc_store(id="simulated-run", data=Dict()),           # the one that gets popped
                dcc_store(id="sim-load", data = 0 ),  # Initially off
                html_br(),
                html_div(id="simloader-msg", children=[]),
                html_br(),
                html_h3("Run animation"),
                html_button("Star/Stop", id="animation-button", n_clicks=0),
                dcc_store(id="animation-on", data=false)  # Initially off
                ])
            ]),
        html_div(id="blockbelow", 
            children = 
                [ 
                dcc_graph(id = "graph",  figure=PlotlyJS.plot(), animate = true,
                    style = Dict(
                    :position => "relative",    
                    :left => "30px",
                    :top => "650px" ))
                ]),
                dcc_checklist(id ="VariablesPlot",
                    options =[
                        Dict("label" => "Assimilates", "value" => "A"),
                        Dict("label" => "T1", "value" => "T1"),
                        Dict("label" => "T2", "value" => "T2"),
                        Dict("label" => "CO2", "value" => "CO2"),
                        Dict("label" => "VP", "value" => "VP"),
                        Dict("label" => "RH", "value" => "RH"),
                        Dict("label" => "Tp", "value" => "Tp"),
                        Dict("label" => "Tout", "value" => "Tout"),
                        Dict("label" => "VPout", "value" => "VPout"),
                        Dict("label" => "RHout", "value" => "RHout"),
                        Dict("label" => "PAR", "value" => "PAR"),
                        Dict("label" => "SW", "value" => "SW"),
                        Dict("label" => "Wind", "value" => "Wind"),
                        Dict("label" => "Tsky", "value" => "Tsky"),
                        Dict("label" => "Tsoil", "value" => "Tsoil")
                    ], style = Dict(
                    :color=> "white", 
                    :position => "relative",    
                    :left => "1100px",
                    :top => "300px"),
                    value=["T1", "T2"]
                )
                
end


# Toggle simulation load variable
callback!(
    app,
    Output("sim-load", "data"),
    Input("sim-button", "n_clicks"),
    State("sim-load", "data")
    ) do clicks, current_state
        # Toggle boolean value
        new_state = !(current_state==1)
        return new_state
end

# Toggle animation-on variable
callback!(
    app,
    Output("animation-on", "data"),
    Input("animation-button", "n_clicks"),
    State("animation-on", "data")
    ) do clicks, current_state
        # Toggle boolean value
        new_state = !(current_state==1)
        return new_state
end

# Simulation enable/disable with start/stop button
callback!(
    app,
    Output("interval", "disabled"),
    Input("animation-on", "data")
    ) do is_on
        return is_on == 1  # If animation is ON, interval should be ENABLED (disabled = false)
end


# Callback simulation load
callback!(
    app,
    Output("simloader-msg","children"),
    Output("simulated-run", "data"),
    Input("sim-load", "data"),
    State("dias_sim", "value"),
    State("date-picker","date")
    ) do simtoogle, noDias, dia
        noDias = noDias === nothing ? 1 : noDias  # Default to 1 day if no value
        if simtoogle === 0
            full_data = Dict()
            return  "No simulation loaded", full_data
        else
            full_data = HarverstSim.harvestSimulation!(Dates.DateTime(dia), noDias, deepcopy(data_history))
            safe_data = clean_dict(full_data)
            return  "Simulation loaded", deepcopy(safe_data)
        end
end

callback!(# Overlay
    app,
    Output("overlay", "children"),
    Input("interval", "n_intervals"),
    State("simulated-run", "data"),    # <-- read from store
    State("sim-load", "data"),
    ) do n, data, simtoogle

    if !isempty(data) && (simtoogle !== 1)
        dd = deep_collect(data)
    
        if length(dd["T1"]) > n
            T1 =   dd["T1"][n]
            T2 =   dd["T2"][n]
            CO2 =   dd["CO2"][n]
            RH =   dd["RH"][n]
            VP =   dd["VP"][n]
            PAR =   dd["PAR"][n]
            Tp =   dd["Tp"][n]
            SW =   dd["SW"][n]
            Tout =   dd["Tout"][n]
            Tsky =   dd["Tsky"][n]
            Tsoil =   dd["Tsoil"][n]
            Wind =   dd["Wind"][n]
            RHout =   dd["RHout"][n]
            VPout =   dd["VPout"][n]
        else
            T1 = T2 = CO2 = RH = VP = PAR = Tp = SW = Tout = Tsky = Tsoil = Wind = RHout = VPout = 0.0
        end

        overlay = [
            create_text_block(T1, "T1", "rgb(133,178,133)", "165px", "460px", "C"),
            create_text_block(T2,"T2", "cornflowerblue", "320px", "280px", "C"),
            create_text_block(CO2,"CO2", "cornflowerblue", "790px", "360px", "mg/m^3"),
            create_text_block(RH,"RH", "cornflowerblue", "600px", "280px", "%"), 
            create_text_block(VP,"VP", "cornflowerblue", "740px", "280px", "Pa"), 
            create_text_block(PAR,"PAR", "orange", "480px", "365px", "W/m^2"),
            create_text_block(Tp,"Tp", "rgb(238,128,124)", "750px", "592px", "C"),
            create_text_block(SW,"SW", "gold", "800px", "190px", "W/m^2"),
            create_text_block(Tout,"Tout", "skyblue", "175px", "190px", "C"),
            create_text_block(Tsky,"Tsky", "lightskyblue ", "175px", "15px", "C"),
            create_text_block(Tsoil,"Tsoil", "lightgray", "175px", "592px", "C"),
            create_text_block(Wind,"Wind", "skyblue", "100px", "100px", "m/s"),
            create_text_block(RHout,"RHout", "skyblue", "620px", "190px", "%"), 
            create_text_block(VPout,"VPout", "skyblue", "620px", "160px", "Pa") 
        ]

        return overlay
    else 
        return []
    end
end


callback!(# graph
    app,
    Output("graph", "figure"),
    Input("sim-load","data"),
    Input("VariablesPlot","value"),
    State("dias_sim", "value"),
    State("simulated-run", "data")
    #State("sim-load", "data")
 ) do toggle, variables, dias, data 

    dd2 = deep_collect(data)
    traces = AbstractTrace[] 
    if (toggle !== true) && !isempty(data) 
        for var in variables
            if haskey(dd2, var) && length(dd2[var]) > 1
                y = dd2[var]
                x = 1:length(y)
                push!(traces, scatter(
                    x = x,
                    y = y,
                    mode = "lines",
                    name = var#,
                    #line = attr(color = "#7FDBFF")  # You can use a color map here too
                ))
            end
        end
    end
    layout = Layout(
        title = "Greenhouse Variables",
        paper_bgcolor =  "rgb(10,10,10,0.8)",
        plot_bgcolor = "rgb(10,10,10,0.8)",
        font = attr(color = "#7FDBFF"),
        xaxis = attr(
            color = "#7FDBFF",
            gridcolor = "gray",
            zerolinecolor = "gray"
        ),
        yaxis = attr(
            color = "#7FDBFF",
            gridcolor = "gray",
            zerolinecolor = "gray"
        )
    )


    fig = PlotlyJS.plot(traces, layout)

    return  fig

end



run_server(app, "0.0.0.0", 8050, debug = true)

