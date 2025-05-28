import React, { useState } from 'react';
import Papa from 'papaparse';
import {
  ChakraProvider,
  Box,
  VStack,
  Heading,
  Text,
  Input,
  Button,
  Table,
  Thead,
  Tbody,
  Tr,
  Th,
  Td,
  Alert,
  AlertIcon,
  FormControl,
  FormLabel,
  NumberInput,
  NumberInputField,
  Stat,
  StatLabel,
  StatNumber,
  StatGroup,
} from '@chakra-ui/react';
import { Line } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';
import { optimizePortfolio } from './utils/optimizer';

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
);

function App() {
  const [prices, setPrices] = useState(null);
  const [weights, setWeights] = useState(null);
  const [error, setError] = useState(null);
  const [targetReturn, setTargetReturn] = useState(0.0005);
  const [chartData, setChartData] = useState(null);
  const [portfolioMetrics, setPortfolioMetrics] = useState(null);

  const handleFileUpload = (event) => {
    const file = event.target.files[0];
    if (file) {
      Papa.parse(file, {
        complete: (results) => {
          try {
            const data = results.data;
            if (data.length < 2) {
              throw new Error('CSV file must contain at least two rows');
            }

            const headers = data[0].filter(h => h.trim() !== '');
            if (headers.length < 2) {
              throw new Error('CSV must contain at least two assets');
            }

            const priceData = data.slice(1)
              .map(row => row.slice(0, headers.length)
                .map(val => parseFloat(val)))
              .filter(row => row.every(val => !isNaN(val)));

            if (priceData.length < 2) {
              throw new Error('CSV must contain at least two valid price rows');
            }

            setPrices({ headers, data: priceData });
            setError(null);

            const chartData = {
              labels: Array.from({ length: priceData.length }, (_, i) => i + 1),
              datasets: headers.map((header, idx) => ({
                label: header,
                data: priceData.map(row => row[idx]),
                borderColor: getRandomColor(),
                tension: 0.1
              }))
            };
            setChartData(chartData);
          } catch (err) {
            setError(err.message);
            setPrices(null);
            setWeights(null);
            setChartData(null);
            setPortfolioMetrics(null);
          }
        }
      });
    }
  };

  const calculateWeights = () => {
    if (!prices) return;

    try {
      const result = optimizePortfolio(prices.data, targetReturn);
      
      setWeights(result.weights.map((w, i) => ({
        asset: prices.headers[i],
        weight: w
      })));
      
      setPortfolioMetrics({
        expectedReturn: result.metrics.expectedReturn,
        volatility: result.metrics.volatility
      });

      setError(null);
    } catch (err) {
      setError('Error calculating weights: ' + err.message);
      setWeights(null);
      setPortfolioMetrics(null);
    }
  };

  const getRandomColor = () => {
    const letters = '0123456789ABCDEF';
    let color = '#';
    for (let i = 0; i < 6; i++) {
      color += letters[Math.floor(Math.random() * 16)];
    }
    return color;
  };

  return (
    <ChakraProvider>
      <Box p={8}>
        <VStack spacing={6} align="stretch">
          <Heading>Portfolio Optimizer</Heading>
          
          <Text>
            Upload a CSV file with asset prices. The first row should contain asset names,
            and subsequent rows should contain prices for each asset.
          </Text>

          <FormControl>
            <FormLabel>Upload Price Data (CSV)</FormLabel>
            <Input
              type="file"
              accept=".csv"
              onChange={handleFileUpload}
              padding={1}
            />
          </FormControl>

          <FormControl>
            <FormLabel>Target Return (per period)</FormLabel>
            <NumberInput
              value={targetReturn}
              onChange={(_, val) => setTargetReturn(val)}
              step={0.0001}
              precision={4}
            >
              <NumberInputField />
            </NumberInput>
          </FormControl>

          {error && (
            <Alert status="error">
              <AlertIcon />
              {error}
            </Alert>
          )}

          {prices && (
            <Button colorScheme="blue" onClick={calculateWeights}>
              Calculate Optimal Weights
            </Button>
          )}

          {portfolioMetrics && (
            <StatGroup>
              <Stat>
                <StatLabel>Expected Return (per period)</StatLabel>
                <StatNumber>{(portfolioMetrics.expectedReturn * 100).toFixed(2)}%</StatNumber>
              </Stat>
              <Stat>
                <StatLabel>Portfolio Volatility (per period)</StatLabel>
                <StatNumber>{(portfolioMetrics.volatility * 100).toFixed(2)}%</StatNumber>
              </Stat>
            </StatGroup>
          )}

          {chartData && (
            <Box h="400px">
              <Line
                data={chartData}
                options={{
                  responsive: true,
                  maintainAspectRatio: false,
                  plugins: {
                    title: {
                      display: true,
                      text: 'Asset Prices Over Time'
                    }
                  }
                }}
              />
            </Box>
          )}

          {weights && (
            <Table variant="simple">
              <Thead>
                <Tr>
                  <Th>Asset</Th>
                  <Th isNumeric>Weight</Th>
                  <Th isNumeric>Percentage</Th>
                </Tr>
              </Thead>
              <Tbody>
                {weights.map(({ asset, weight }) => (
                  <Tr key={asset}>
                    <Td>{asset}</Td>
                    <Td isNumeric>{weight.toFixed(4)}</Td>
                    <Td isNumeric>{(weight * 100).toFixed(2)}%</Td>
                  </Tr>
                ))}
              </Tbody>
            </Table>
          )}
        </VStack>
      </Box>
    </ChakraProvider>
  );
}

export default App;