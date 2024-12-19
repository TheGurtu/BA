using Pkg
Pkg.activate(".")
using Agents, Graphs, Random, GraphPlot,Compose, Cairo,Colors, Plots, Distributions, DataFrames, StatsPlots, BenchmarkTools,CSV


@agent struct PersonAgent(GraphAgent)
    S::Bool
    I::Bool
    R::Bool
	infected_days::Int8
	risk_perception::Float32
	rp_probability::Float32
end

function person_step!(agent, model)

	# Wenn der Agent susceptible ist
	if agent.S == true
		infected_neighbors = 0
		nb = 0
		for neighbor_id in neighbors(model.net,agent.id)
			neighbor = model[neighbor_id]
			nb = nb + 1
			if neighbor.I == true
				infected_neighbors = infected_neighbors +1
			end
		end
		#Berechnet Wahrscheinlichkeit fuer Schutzverhalten
		rp = exp(-1*agent.risk_perception*(infected_neighbors/nb))
		probability = (1-(1-(rp*model.pb))^infected_neighbors)*100
		agent.rp_probability = probability
		random = rand(Uniform(1,100))
		if random <= probability
			#Wechsel zum Zustand infected wenn die WK es hergibt
			agent.S = false
			agent.I = true
			model.membership[agent.id] = 2
			return
		end
		
	end
	
	#Wenn der Agent infiziert ist reduziere die infizierten Tage und wechsel Zustand
	if agent.I == true
		if agent.infected_days == 0
			agent.I = false
			agent.R = true
			model.membership[agent.id] = 3
		end
		if agent.infected_days != 0
			agent.infected_days = agent.infected_days -1
		end

		return
	end

	if agent.R == true
		return
	end

	#Wenn Agent susceptible ist naehere die Meinung an die Nachbarn an
	neighbour_with_near_opinion = 0
	opinion_sum = 0
	for neighbor_id in neighbors(model.net, agent.id)
		neighbor = model[neighbor_id]
		if abs(neighbor.risk_perception-agent.risk_perception) <= model.threshold && neighbor.S
			neighbour_with_near_opinion = neighbour_with_near_opinion +1
			opinion_sum = opinion_sum + neighbor.risk_perception - agent.risk_perception

			
		end
		
	end
	
	
	if neighbour_with_near_opinion != 0
		change_rate = opinion_sum/neighbour_with_near_opinion
		agent.risk_perception = agent.risk_perception + (model.mu*change_rate)
		
	end
	


end


function initialize(pb, mu, num_agents, max_k, risk_dist, threshold)
	
	rng = Random.Xoshiro(rand(1:200))
    net = watts_strogatz(num_agents, max_k, 0.1)
    loc_x, loc_y = spring_layout(net)
    space = GraphSpace(net)

    properties = Dict(:pb => pb, :mu => mu, :net => net, :loc_x => loc_x, :loc_y => loc_y, :membership => ones(Int8, num_agents), :num_agents => num_agents, :threshold => threshold)
	model = StandardABM(PersonAgent, space; properties, agent_step! = person_step! , rng, scheduler= Schedulers.Randomly())
	

	for n in 1:num_agents
		agent = PersonAgent(n, n, true, false, false,5,risk_dist[n],pb) 
		if n == 1
			agent = PersonAgent(n, n, false, true, false,5,risk_dist[n],pb)
			model.membership[n] = 2
		end
		add_agent_single!(agent, model)
	end
	return model
end 

function generate_array(num_agents::Int, method::Symbol, min_rp::Int, max_rp::Int, mu::Float64, sigma::Float64)
    if method == :random
        # Zufällige Verteilung (gleichverteilte Werte zwischen 1 und 10)
        return rand(min_rp:max_rp, num_agents)
    elseif method == :normal
		mu = 5.5
		sigma = 1.2
		dist = Normal(mu, sigma)
		risk_dist = rand(dist, num_agents)
        # Normalverteilte Werte (zwischen 1 und 10 begrenzt)
        return rand(dist, num_agents)
    
    elseif method == :bi
        # Gleichmäßig verteilte Werte zwischen 1 und 10 (nur 1 oder 10)
        return [mod(i,2)==1 ? min_rp : max_rp for i in 1:num_agents]
    
    else
        error("Unbekannte Methode: $method")
    end
end

function generate_array(n::Int, const_value::Float64)
    # Array mit konstantem Wert für alle Stellen
    return fill(const_value, n)
end

print("\n")

#Simualtionsinterface fuer die weitere Simulationen
function simulate(num_simualations, max_steps,risk_dist,prefix;num_agents=1000, max_degree=8, pb=0.25, mu=0.5, threshold=10, constant=false, const_val=0, min_rp=1, max_rp=10, ew=5.0,variance=1.2 )
	
	#functions for information extraction
	extract_inf(model) = return count([agent.I for agent in allagents(model)])
	extract_sus(model) = return count([agent.S for agent in allagents(model)])
	extract_re(model) = return count([agent.R for agent in allagents(model)])
	extract_numAgents(model) = return model.num_agents
	extract_mu(model) = return model.mu
	extract_thr(model) = return model.threshold
	#extract_avg_rp(model) return mean([agent.risk_perception for agent in allagents(model)end])
	function extract_avg_rp(model)
		risk_vec = []
		for agent in allagents(model)
			push!(risk_vec, agent.risk_perception)
		end
		return mean(risk_vec)
	end


	model_vector = []
	for i in 1:num_simualations
		if constant == false
			push!(model_vector,  initialize(pb,mu, num_agents, max_degree, generate_array(num_agents, risk_dist, min_rp, max_rp, ew, variance), threshold))
		else
			push!(model_vector,  initialize(pb,mu, num_agents, max_degree, generate_array(num_agents, const_val), threshold))
		end
	end

	complete_model = DataFrame(time=[],extract_inf=[],extract_sus=[], extract_re=[], extract_avg_rp=[],extract_numAgents=[], extract_mu=[], extract_thr=[], ensemble=[], peak_time=[], duration = [])
	complete_agents = DataFrame(time=[],id=[], risk_perception=[],S=[], I=[], R=[],infected_days=[])
	for i in 1:num_simualations
		df_agent, df_model = Agents.run!(model_vector[i], max_steps, adata=[:risk_perception, :S, :I, :R, :infected_days], mdata=[extract_inf, extract_sus, extract_re, extract_avg_rp, extract_numAgents, extract_mu, extract_thr])
		insertcols!(df_model,:ensemble => i)
		
		max = maximum(df_model[:,:extract_inf])
		
		max_time = filter(:extract_inf => inf->inf==max, df_model)
	
		insertcols!(df_model,:peak_time => max_time.time[1])
		duration = filter(:extract_inf => d->d==0,df_model)
		insertcols!(df_model,:duration => duration.time[1])

		complete_model = vcat(complete_model, df_model)
		complete_agents = vcat(complete_agents, df_agent)
	end


	
	insertcols!(complete_model, :steps => max_steps)
	insertcols!(complete_model, :ident => prefix)
	insertcols!(complete_agents, :ident => prefix)

	return complete_model, complete_agents

end



function scenario_1()
	prefix="influence_of_rp"
	num_simualations = 50
	steps = 200

	df_model, df_agent = simulate(num_simualations, steps,:constant, prefix*"_no_rp"; constant=true, const_val=0.0)
	normal_dist = simulate(num_simualations,steps,:normal, prefix*"_normal_dist_ew5")
	df_model = vcat(df_model, normal_dist[1])
	df_agent = vcat(df_agent, normal_dist[2])
	return df_model, df_agent

end


function scenario_2()
	prefix="compare_different_constant_rp"
	num_simualations = 50
	steps = 200
	
	df_model, df_agent = simulate(num_simualations, steps,:constant, prefix*"_constant_rp_1"; constant=true, const_val=1.0)

	rp_2 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_2"; constant=true, const_val=2.0)
	df_model = vcat(df_model, rp_2[1])
	df_agent = vcat(df_agent, rp_2[2])

	rp_3 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_3"; constant=true, const_val=3.0)
	df_model = vcat(df_model, rp_3[1])
	df_agent = vcat(df_agent, rp_3[2])

	rp_4 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_4"; constant=true, const_val=4.0)
	df_model = vcat(df_model, rp_4[1])
	df_agent = vcat(df_agent, rp_4[2])

	CSV.write("scenario_2_1_model.csv", df_model)
	CSV.write("scenario_2_1_agent.csv", df_agent)

	df_model, df_agent = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_5"; constant=true, const_val=5.0)


	rp_6 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_6"; constant=true, const_val=6.0)
	df_model = vcat(df_model, rp_6[1])
	df_agent = vcat(df_agent, rp_6[2])

	rp_7 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_7"; constant=true, const_val=7.0)
	df_model = vcat(df_model, rp_7[1])
	df_agent = vcat(df_agent, rp_7[2])

	rp_8 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_8"; constant=true, const_val=8.0)
	df_model = vcat(df_model, rp_8[1])
	df_agent = vcat(df_agent, rp_8[2])
	CSV.write("scenario_2_2_model.csv", df_model)
	CSV.write("scenario_2_2_agent.csv", df_agent)

	df_model, df_agent = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_9"; constant=true, const_val=9.0)


	rp_10 = simulate(num_simualations,steps,:constant, prefix*"_constant_rp_10"; constant=true, const_val=10.0)
	df_model = vcat(df_model, rp_10[1])
	df_agent = vcat(df_agent, rp_10[2])
	CSV.write("scenario_2_3_model.csv", df_model)
	CSV.write("scenario_2_3_agent.csv", df_agent)


end

function scenario_3()
	prefix="compare_different_dist_rp_no_threshold"
	num_simualations = 50
	steps = 200

	df_model, df_agent = simulate(num_simualations,steps,:bi, prefix*"_bi_dist_rp")
	rand_dist = simulate(num_simualations,steps,:random, prefix*"_random_dist")
	df_model = vcat(df_model,rand_dist[1])
	df_agent = vcat(df_agent, rand_dist[2])

	return df_model, df_agent

end


function scenario_5()
	prefix = "compare_different_dist_diff_mu"
	num_simualations = 50
	steps = 200

	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_mu_01"; mu=0.1)
	bi_dist = simulate(num_simualations, steps, :bi, prefix*"_bi_mu_01"; mu=0.1)
	df_model = vcat(df_model, bi_dist[1])
	df_agent = vcat(df_agent, bi_dist[2])
	rnd_dist = simulate(num_simualations, steps, :random, prefix*"_random_mu_01"; mu=0.1)
	df_model = vcat(df_model, rnd_dist[1])
	df_agent = vcat(df_agent, rnd_dist[2])
	CSV.write("scenario_5_1_model.csv", df_model)
	CSV.write("scenario_5_1_agent.csv", df_agent)

	
	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_mu_03"; mu=0.3)

	bi_dist = simulate(num_simualations, steps, :bi, prefix*"_bi_mu_03"; mu=0.3)
	df_model = vcat(df_model, bi_dist[1])
	df_agent = vcat(df_agent, bi_dist[2])
	rnd_dist = simulate(num_simualations, steps, :random, prefix*"_random_mu_03"; mu=0.3)
	df_model = vcat(df_model, rnd_dist[1])
	df_agent = vcat(df_agent, rnd_dist[2])
	CSV.write("scenario_5_2_model.csv", df_model)
	CSV.write("scenario_5_2_agent.csv", df_agent)

	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_mu_07"; mu=0.7)

	bi_dist = simulate(num_simualations, steps, :bi, prefix*"_bi_mu_07"; mu=0.7)
	df_model = vcat(df_model, bi_dist[1])
	df_agent = vcat(df_agent, bi_dist[2])
	rnd_dist = simulate(num_simualations, steps, :random, prefix*"_random_mu_07"; mu=0.7)
	df_model = vcat(df_model, rnd_dist[1])
	df_agent = vcat(df_agent, rnd_dist[2])
	CSV.write("scenario_5_3_model.csv", df_model)
	CSV.write("scenario_5_3_agent.csv", df_agent)

	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_mu_1"; mu=1.0)

	bi_dist = simulate(num_simualations, steps, :bi, prefix*"_bi_mu_1"; mu=1.0)
	df_model = vcat(df_model, bi_dist[1])
	df_agent = vcat(df_agent, bi_dist[2])
	rnd_dist = simulate(num_simualations, steps, :random, prefix*"_random_mu_1"; mu=1.0)
	df_model = vcat(df_model, rnd_dist[1])
	df_agent = vcat(df_agent, rnd_dist[2])
	CSV.write("scenario_5_4_model.csv", df_model)
	CSV.write("scenario_5_4_agent.csv", df_agent)


end




function scenario_4()
	prefix="compare_different_dist_rp_with_threshold"
	num_simualations = 50
	steps = 200


	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_distributed_risk_threshold_7";threshold=7)

	rand_dist = simulate(num_simualations,steps,:random, prefix*"_random_dist_threshold_7";threshold=7)
	df_model = vcat(df_model,rand_dist[1])
	df_agent = vcat(df_agent, rand_dist[2])
	

	normal_dist = simulate(num_simualations, steps, :normal, prefix*"_normal_distributed_risk_threshold_5";threshold=5)
	df_model = vcat(df_model, normal_dist[1])
	df_agent = vcat(df_agent, normal_dist[2])

	rand_dist = simulate(num_simualations,steps,:random, prefix*"_random_dist_threshold_5";threshold=5)
	df_model = vcat(df_model,rand_dist[1])
	df_agent = vcat(df_agent, rand_dist[2])
	CSV.write("scenario_4_1_model.csv", df_model)
	CSV.write("scenario_4_1_agent.csv", df_agent)

	df_model, df_agent = simulate(num_simualations, steps, :normal, prefix*"_normal_distributed_risk_threshold_2";threshold=2)

	rand_dist = simulate(num_simualations,steps,:random, prefix*"_random_dist_threshold_2";threshold=2)
	df_model = vcat(df_model,rand_dist[1])
	df_agent = vcat(df_agent, rand_dist[2])
	CSV.write("scenario_4_2_model.csv", df_model)
	CSV.write("scenario_4_2_agent.csv", df_agent)

	

end

function scenario_6()
	prefix="const_rp_5.5"
	num_simualations = 50
	steps = 200

	df_model, df_agent = simulate(num_simualations,steps,:constant, prefix; constant=true, const_val=5.5)
	
	#insertcols!(df_model, :scenario => prefix)
	return df_model, df_agent

end



function scenario_7()
	prefix="different_dist_no_influence_"
	num_simualations = 50
	steps = 200

	df_model, df_agent = simulate(num_simualations,steps,:normal, prefix*"_normal_dist_ew5"; threshold=0, mu=0.0)
	bi_dist = simulate(num_simualations,steps,:bi, prefix*"_bi_dist_rp"; threshold=0, mu=0.0)
	df_model = vcat(df_model,bi_dist[1])
	df_agent = vcat(df_agent, bi_dist[2])
	rand_dist = simulate(num_simualations,steps,:random, prefix*"_random_dist"; threshold=0, mu=0.0)
	df_model = vcat(df_model,rand_dist[1])
	df_agent = vcat(df_agent, rand_dist[2])
	CSV.write("scenario_7_model.csv", df_model)
	CSV.write("scenario_7_agent.csv", df_agent)
	

end


df_model, df_agent = scenario_1()
CSV.write("scenario_1_model.csv", df_model)
CSV.write("scenario_1_agent.csv", df_agent)
scenario_2()
df_model, df_agent = scenario_3()
CSV.write("scenario_3_model.csv", df_model)
CSV.write("scenario_3_agent.csv", df_agent)
scenario_4()
scenario_5()
df_model, df_agent = scenario_6()
CSV.write("scenario_6_model.csv", df_model)
CSV.write("scenario_6_agent.csv", df_agent)
scenario_7()



